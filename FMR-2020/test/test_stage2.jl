using Test
using FMR2020

@testset "Stage 2 endogenous nonlinear transition blocks" begin
    p=default_parameters(n_coarse=61,n_fine=241,
                         egm_tol=1e-9,dist_tol=1e-9)
    ref=baseline_reference(p)
    chain=rouwenhorst(p)
    grids=make_asset_grids(p,ref.K)
    B=ref.nfa_y*ref.Y
    aggregate=FMR2020.aggregate_prices(p,ref,ref.r,B)
    pr=HouseholdPrices(ref.r,aggregate.wage,aggregate.terms,
                       aggregate.L,aggregate.profits)
    hh=solve_household(p,grids,chain,pr,0.9832)
    mu=stationary_distribution(p,grids,chain,hh).distribution
    st=stationary_candidate(p,ref,chain,grids,0.9832,ref.r,B)
    Tperiods=8
    inp=TransitionInputs(p,ref,chain,grids,st,st,
                         fill(B,Tperiods),ones(Float64,Tperiods),
                         17.0,1.0,B)
    x=vcat(fill(st.K,Tperiods-1),fill(st.L,Tperiods))
    @test length(x)==2*Tperiods-1
    @testset "Credit expansion path and time indexing" begin
        boom=credit_expansion_path(p,ref.Y;periods=12)
        @test length(boom)==12
        @test isapprox(boom[1],-0.074*ref.Y;rtol=1e-12)
        @test isapprox(boom[2]-boom[1],-0.074*ref.Y*0.98;rtol=1e-12)
        @test all(diff(boom) .< 0)
        @test isapprox(boom[end],-0.074*ref.Y*(1-0.98^12)/(1-0.98);rtol=1e-12)
    end
    @testset "Prices and period bounds" begin
        workspace=TransitionWorkspace(inp)
        FMR2020.transition_prices!(workspace,inp,x)
        @test maximum(abs.(workspace.q .- 1)) < 1e-10
        @test maximum(abs.(workspace.K .- st.K)) < 1e-10
        @test maximum(abs.(workspace.L .- st.L)) < 1e-10
        @test minimum(workspace.terms) > 0
        @test all(isfinite,workspace.rate)
        @test_throws DimensionMismatch FMR2020.transition_prices!(
            workspace,inp,x[1:end-1])
    end
    @testset "Backward household EGM and forward distribution" begin
        trial=evaluate_transition!(inp,x)
        w=trial.workspace
        @test length(trial.residual)==15
        @test size(w.coarse_consumption,3)==8
        @test size(w.fine_savings,3)==8
        @test size(w.distribution,3)==8
        @test all(isfinite,trial.residual)
        @test minimum(w.fine_consumption)>0
        @test minimum(w.fine_savings) >= grids.borrowing_limit
        @test maximum(w.distribution_error) < 1e-8
        @test isapprox(sum(view(w.distribution,:,:,size(w.distribution,3))),1;atol=1e-8)
        @test all(isfinite,w.terms)
        @test maximum(abs.(w.q .- 1)) < 1e-10
        result=transition_stationary_check(inp)
        @test isfinite(result.max_scaled_residual)
        @test result.max_scaled_labor >= 0
        # The author's code omits the power 1/varphi in labor supply.
        # This is an expressly documented source inconsistency, not
        # a reason to silently change the replicated equation.
    end
    @testset "Off-grid saving at unoccupied states is not clipped" begin
        grid=Float64[0.0,1.0]
        P=Float64[1.0 0.0;0.0 1.0]
        policy=Float64[0.5 0.5;1.1 1.2]
        mass=zeros(2,2)
        mass[1,1]=1.0
        next_mass=similar(mass)
        distribution_step!(next_mass,mass,policy,grid,P)
        @test isapprox(sum(next_mass),1.0;atol=1e-12)
        @test isapprox(next_mass[1,1],0.5;atol=1e-12)
        @test isapprox(next_mass[2,1],0.5;atol=1e-12)
        mass[1,1]=0.0
        mass[2,1]=1.0
        @test_throws DomainError distribution_step!(next_mass,mass,policy,grid,P)
    end

    @testset "Independent interior goods vs terminal boundary discrepancy" begin
        trial=evaluate_transition!(inp,x)
        w=trial.workspace
        assumed=FMR2020.TransitionSolveResult(inp,w,copy(x),trial.residual,
                                                0,1,maximum(abs,trial.residual),true)
        check=transition_diagnostics(assumed;compute_euler=false)
        @test check.max_goods_market_period in 0:(Tperiods-1)
        @test isfinite(check.max_goods_market_level)
        @test isfinite(check.terminal_goods_gap_unconstrained)
        @test isfinite(check.terminal_wealth_gap_unconstrained)
        @test check.max_goods_market_level >= 0
        # Only the last-period household aggregate is changed. The
        # author model does not impose a T-th goods-market condition.
        last_home=w.aggregate_home[end]
        w.aggregate_home[end]=last_home+0.25
        modified=transition_diagnostics(assumed;compute_euler=false)
        @test isapprox(modified.max_goods_market_level,
                       check.max_goods_market_level;atol=1e-12)
        @test isapprox(modified.terminal_goods_gap_unconstrained,
                       check.terminal_goods_gap_unconstrained-0.25;atol=1e-12)
        w.aggregate_home[end]=last_home
    end

    @testset "Export guard and solver input validation" begin
        @test_throws DimensionMismatch evaluate_transition!(inp,x[1:end-1])
    end
end
