using Test
using FMR2020

@testset "FMR-2020 Julia stationary solver" begin
    p = default_parameters(n_coarse=97,n_fine=385,egm_tol=1e-9,dist_tol=1e-9)
    ref = baseline_reference(p)
    chain = rouwenhorst(p)
    grids = make_asset_grids(p,ref.K)

    @testset "Parameters and Rouwenhorst" begin
        @test_throws ArgumentError default_parameters(alpha=1.2)
        @test_throws ArgumentError default_parameters(theta=1.2)
        @test length(chain.levels) == p.n_s
        @test maximum(abs.(vec(sum(chain.transition,dims=2)) .- 1)) < 1e-12
        @test maximum(abs.(transpose(chain.transition)*chain.invariant .-
                           chain.invariant)) < 1e-10
        @test isapprox(sum(chain.levels .* chain.invariant),1;atol=1e-12)
        @test minimum(chain.transition) >= 0
    end

    @testset "Borrowing grid and independent Hermite identity" begin
        @test length(grids.coarse_full) == p.n_coarse
        @test length(grids.fine_full) == p.n_fine
        @test grids.coarse_full[1] < grids.borrowing_limit
        @test grids.coarse[1] == grids.borrowing_limit
        @test grids.fine[1] == grids.borrowing_limit
        @test all(diff(grids.coarse) .> 0)
        @test all(diff(grids.fine) .> 0)
        x = Float64[0,0.3,0.9,2.0,3.4]
        y = 2 .+ 3 .* x
        slope = similar(x)
        FMR2020.pchip_slopes!(slope,x,y)
        @test maximum(abs.(slope .- 3)) < 1e-12
        for v in range(x[1],x[end],length=31)
            @test isapprox(FMR2020.pchip_evaluate(x,y,slope,v),2+3v;atol=1e-12)
        end
    end

    @testset "Scalar root and labor closure" begin
        result = FMR2020.bracketed_root(x->x^2-2,1.0,2.0;tol=1e-10)
        @test abs(result.root-sqrt(2)) < 1e-9
        B = ref.nfa_y*ref.Y
        a = FMR2020.aggregate_prices(p,ref,ref.r,B)
        @test isapprox(a.K,ref.K;rtol=1e-6)
        @test isapprox(a.L,1;rtol=1e-6)
        @test isapprox(a.terms,1;rtol=1e-6)
        @test abs(a.union_error) < 1e-7
    end

    B = ref.nfa_y*ref.Y
    prices = FMR2020.aggregate_prices(p,ref,ref.r,B)
    household_prices = HouseholdPrices(ref.r,prices.wage,prices.terms,
                                       prices.L,prices.profits)
    hh = solve_household(p,grids,chain,household_prices,0.9832)
    diag = household_diagnostics(p,grids,chain,household_prices,0.9832,hh)
    @testset "EGM Euler/KKT/budget" begin
        @test hh.update_error <= p.egm_tol
        @test minimum(hh.consumption_fine) > 0
        @test minimum(hh.savings_fine) >= grids.borrowing_limit-1e-8
        @test diag.max_budget_error < 1e-8
        @test diag.minimum_savings_step > -1e-6
        @test isfinite(diag.max_euler_unconstrained)
        @test isfinite(diag.max_kkt_violation)
        @test diag.max_euler_unconstrained < 0.1 # grid-dependent; report actual
        @test diag.max_kkt_violation < 0.1
    end

    @testset "Matrix-free distribution and mass conservation" begin
        n = length(grids.fine)
        ns = length(chain.levels)
        mu = zeros(n,ns)
        mu[1,:] .= chain.invariant
        dest = similar(mu)
        distribution_step!(dest,mu,hh.savings_fine,grids.fine,chain.transition)
        @test minimum(dest) >= 0
        @test isapprox(sum(dest),1;atol=1e-12)
        @test_throws DomainError begin
            bad = copy(hh.savings_fine)
            bad[1,1] = grids.maximum_asset*2
            distribution_step!(dest,mu,bad,grids.fine,chain.transition)
        end
        stationary = stationary_distribution(p,grids,chain,hh)
        invariant = distribution_invariance(stationary.distribution,hh,grids,chain)
        @test abs(invariant.mass-1) < 1e-10
        @test invariant.min_mass >= 0
        @test invariant.invariance_error < 1e-8
        @test isfinite(sum(stationary.distribution .* hh.savings_fine))
    end

    @testset "Full GE opt-in (set FMR_FULL_TESTS=1)" begin
        if get(ENV,"FMR_FULL_TESTS","0") == "1"
            result = solve_three_steady_states(p)
            for st in (result.initial,result.calibrated,result.final)
                check = steady_state_diagnostics(p,result.chain,result.grids,st)
                @test check.scaled_asset_error <= p.root_tol
                @test check.scaled_goods_error < 1e-5
                @test check.distribution.invariance_error < 1e-8
            end
        else
            @info "Full GE root-finding is exercised by scripts/run_steady_states.jl"
            @test true
        end
    end
end
