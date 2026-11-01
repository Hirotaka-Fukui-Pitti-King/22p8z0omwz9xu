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

    @testset "Feasible interior interest-rate bracketing" begin
        # Regression: a fixed upper endpoint fails for an unrelated rate even
        # when a root exists at a perfectly feasible interior value.
        function restricted_residual(r::Float64)
            r > 0.014 && throw(DomainError(r,
                "fine-grid policy exceeds distribution grid upper bound"))
            return r-0.011
        end
        bracket = FMR2020.equilibrium_rate_bracket(restricted_residual,
                                                    0.010,0.0,0.017;
                                                    tolerance=1e-12)
        @test bracket.lower <= 0.011 <= bracket.upper
        result = FMR2020.bracketed_root(restricted_residual,
                                        bracket.lower,bracket.upper;
                                        tol=1e-12)
        @test abs(result.root-0.011) < 1e-9
        @test_throws ArgumentError FMR2020.equilibrium_rate_bracket(
            r->r+1.0,0.01,0.0,0.017;maximum_steps=4)
        @test_throws DomainError FMR2020.equilibrium_rate_bracket(
            r->throw(DomainError(r,"unexpected economic error")),0.01,0.0,0.017)
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
        @test diag.max_fine_budget_error < 1e-8
        @test diag.minimum_fine_saving >= grids.borrowing_limit
        @test diag.minimum_savings_step > -1e-6
        # The fine-grid interpolation preserves coarse policy choices at
        # shared knots (full and active grids are nested for these sizes).
        for (i, a) in enumerate(grids.coarse)
            jfine = findfirst(==(a),grids.fine)
            if jfine !== nothing
                @test maximum(abs.(hh.savings_fine[jfine,:] .-
                                    hh.savings_coarse[i,:])) < 1e-8
            end
        end
        @test isfinite(diag.max_euler_unconstrained)
        @test isfinite(diag.max_kkt_violation)
        @test diag.max_euler_unconstrained < 0.1 # grid-dependent; report actual
        @test diag.max_kkt_violation < 0.1
    end

    @testset "Regression: finer grids and impatient beta do not create artificial borrowing" begin
        # The reported failure occurred at this stage, before beta root
        # evaluation: fine-grid savings were about 1e-5 below the limit.
        for (nc,nf) in ((61,241),(193,769))
            pg = default_parameters(n_coarse=nc,n_fine=nf,
                                    egm_tol=1e-9,dist_tol=1e-9)
            grid = make_asset_grids(pg,ref.K)
            sol = solve_household(pg,grid,chain,household_prices,0.980)
            metric = household_diagnostics(pg,grid,chain,household_prices,0.980,sol)
            @test minimum(sol.savings_fine) >= grid.borrowing_limit
            @test metric.max_fine_budget_error < 1e-8
            @test minimum(sol.consumption_fine) > 0
        end
    end


    @testset "Author direct endogenous EGM fine-grid evaluation" begin
        # Regression on the author's numerical algorithm, not on an alternative.
        ng = length(grids.coarse)
        ns = length(chain.levels)
        asset_knots = Vector{Float64}(undef,ng)
        consumption_knots = similar(asset_knots)
        deriv = similar(asset_knots)
        gross = 1+household_prices.r
        price = FMR2020.consumption_cost(household_prices.terms,p)
        for j in 1:ns
            income = FMR2020.labor_income(household_prices,chain.levels[j]) +
                     household_prices.profits
            for i in 1:ng
                expected_mu = sum(chain.transition[j,k] *
                                  hh.consumption_coarse[i,k]^(-p.gamma)
                                  for k in 1:ns)
                consumption_knots[i] = (0.9832*gross*expected_mu)^(-inv(p.gamma))
                asset_knots[i] =
                    (price*consumption_knots[i]+grids.coarse[i]-income)/gross
            end
            FMR2020.pchip_slopes!(deriv,asset_knots,consumption_knots)
            for i in eachindex(grids.fine)
                a=grids.fine[i]
                expected_c = if a <= asset_knots[1]
                    (gross*a+income-grids.borrowing_limit)/price
                else
                    FMR2020.pchip_evaluate(asset_knots,consumption_knots,deriv,a)
                end
                @test isapprox(hh.consumption_fine[i,j],expected_c;atol=1e-10,rtol=1e-10)
            end
        end
        @test minimum(hh.savings_fine) >= grids.borrowing_limit
        @test diag.max_fine_budget_error < 1e-8
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

    @testset "Fine-grid, probability-weighted Euler audit" begin
        st = stationary_candidate(p,ref,chain,grids,0.9832,ref.r,B;
                                  name="accuracy-test")
        result = euler_accuracy_map(p,grids,chain,st)
        s = result.summary
        @test length(result.records) == length(grids.fine)*p.n_s
        @test isapprox(s.total_mass,1;atol=1e-10)
        @test s.n_active_states > 0
        @test s.n_active_states <= s.n_grid_states
        @test s.mass_in_top_ten_percent_asset_grid >= 0
        @test s.probability_weighted_mean >= 0
        @test isfinite(s.max_unconstrained)
        @test isfinite(s.probability_weighted_p99)
        @test result.ranked[1].euler_kkt_violation >=
              result.ranked[end].euler_kkt_violation
        @test maximum(x->x.euler_kkt_violation,result.records) >=
              s.probability_weighted_mean
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
