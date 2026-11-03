using Test
using FMR2020

include(joinpath(@__DIR__,"..","scripts","run_economic_replication.jl"))

@testset "Saved economic market-root resumption is checked, not trusted" begin
    p=default_parameters(n_coarse=61,n_fine=241,
                         egm_tol=1e-9,dist_tol=1e-9)
    steady=solve_three_steady_states(p;verbose=false)
    ctx=stage3_expansion_context(steady)
    n=length(ctx.economy.B)
    @test n==200
    K=fill(ctx.economy.initial.K,n)
    L=fill(ctx.economy.initial.L,n)
    terms=fill(ctx.economy.initial.terms,n)
    mktempdir() do d
        csv=joinpath(d,"synthetic.csv")
        function export_path(;wrong_B=false)
            open(csv,"w") do io
                println(io,
                    "period,B,K,L,Y,I,q,r,w,terms,foreign_consumption,home_consumption,wealth,nominal_exchange_rate")
                for t in 1:n
                    B=ctx.economy.B[t]+((wrong_B && t==2) ? 0.25 : 0.0)
                    println(io,join((t,B,K[t],L[t],1.0,0.1,1.0,0.01,0.5,
                                     terms[t],1.0,1.0,1.0,1.0),","))
                end
            end
        end
        export_path()
        v=read_economic_root_guess(csv,ctx)
        @test length(v)==3n-1
        @test v[1:n-1]==K[2:end]
        @test v[n:2n-1]==L
        @test v[2n:3n-1]==terms
        copy_to=joinpath(d,"economic_synthetic_MARKET_ROOT_UNCERTIFIED.csv")
        cp(csv,copy_to)
        @test baseline_guess_if_present(d,"synthetic",ctx,zeros(3n-1))==v
        @test baseline_guess_if_present(d,"absent",ctx,v)==v
        # The same forward-horizon helper must accept a corrected
        # 200-date market root as its seed, without altering the first
        # 200 dates of B or the aggregate K/L/terms guess.
        extended=FMR2020.stage3_extend_horizon(ctx,204)
        @test extended.economy.B[1:n]==ctx.economy.B
        @test all(extended.economy.B[n+1:end] .== ctx.economy.terminal.B)
        prolonged=FMR2020.stage3_prolong_aggregate_guess(v,ctx,extended)
        @test length(prolonged)==3*204-1
        @test prolonged[1:n-1]==v[1:n-1]
        @test prolonged[204:204+n-1]==v[n:2n-1]
        @test prolonged[2*204:2*204+n-1]==v[2n:3n-1]
        fake=(K=copy(K),L=copy(L),Y=ones(n),I=fill(0.1,n),
              wage=fill(0.5,n),rate=fill(0.01,n),
              terms=copy(terms),aggregate_home=ones(n))
        same=compare_boom_horizons(csv,fake)
        @test same.maximum==0.0
        @test same.first_periods==40
        @test compare_economic_baseline_if_present(
            d,"synthetic",fake).maximum==0.0
        @test compare_economic_baseline_if_present(
            d,"absent",fake)===nothing
        fake.K[5]+=0.03
        different=compare_boom_horizons(csv,fake)
        @test isapprox(different.by_variable["K"],
                       0.03/max(1.0,abs(K[5]));atol=1e-12)
        export_path(wrong_B=true)
        @test_throws ArgumentError read_economic_root_guess(csv,ctx)
        rm(csv)
        @test read_economic_root_guess(csv,ctx)===nothing
    end
    # A four-date steady environment exercises the ACTUAL backward
    # Stage3 household map, the adaptive Euler corrector and the
    # forward transport. No full GE Jacobian or new 599-root calculation.
    Ttiny=4
    base=ctx.economy
    short=TransitionInputs(p,base.reference,base.chain,base.grids,
                           base.initial,base.initial,
                           fill(base.initial.B,Ttiny),ones(Ttiny),
                           base.phi_ac,base.q_initial,base.B_initial)
    shortctx=Stage3Context(short,copy(ctx.mu_pre),ctx.terms_pre,
                           0.0,0.0,true,"4-date EGM accuracy smoke")
    shortx=vcat(fill(base.initial.K,Ttiny-1),
                fill(base.initial.L,Ttiny),
                fill(base.initial.terms,Ttiny))
    work=TransitionWorkspace(short)
    residual=zeros(3Ttiny-1)
    FMR2020.stage3_flexible_full_residual!(
        residual,shortx,shortctx,work)
    polishing=FMR2020.stage3_euler_polish_households!(
        work,shortctx;trigger=1e-5,target=1e-9)
    @test polishing.refined_fine>=0
    @test polishing.refined_coarse>=0
    @test polishing.maximum_euler_violation_on_supported_grid<=1e-5
    FMR2020.stage3_distribution_forward!(work,shortctx)
    @test isapprox(sum(view(work.distribution,:,:,Ttiny)),1.0;atol=1e-9)
    refined_res=similar(residual)
    FMR2020.stage3_flexible_full_market_residual!(
        refined_res,work,shortctx)
    @test all(isfinite,refined_res)
    # This equivalence is the key GE correctness property: the policy
    # correction must be included inside EVERY aggregate residual call.
    work2=TransitionWorkspace(short)
    integrated_res=similar(residual)
    FMR2020.stage3_flexible_full_residual!(
        integrated_res,shortx,shortctx,work2;
        polish_households=true)
    @test integrated_res≈refined_res atol=1e-10 rtol=1e-10
    # NaN at infeasible zero-mass negative-asset nodes is intentional.
    @test isequal(work2.fine_savings,work.fine_savings)
    @test work2.coarse_consumption≈work.coarse_consumption atol=1e-10 rtol=1e-10
    @test work2.aggregate_home≈work.aggregate_home atol=1e-10 rtol=1e-10
    @test all(isfinite,integrated_res)
    # The corrected gross-portfolio map must be selected in the
    # aggregate residual/Jacobian operator, NOT just post-solve.
    # This small 4-period test avoids a heavy new GE root calculation.
    strictwork=TransitionWorkspace(short)
    strictres=similar(residual)
    FMR2020.stage3_flexible_full_residual!(
        strictres,shortx,shortctx,strictwork;
        polish_households=true,strict_portfolio=true)
    @test all(isfinite,strictres)
    @test isapprox(sum(view(strictwork.distribution,:,:,Ttiny)),1.0;
                   atol=1e-9)
    @test_throws ArgumentError FMR2020.stage3_nonlinear_solve(
        shortctx,shortx[1:2Ttiny-1];
        polish_households=true,strict_portfolio=true,
        tolerance=1e-4,maxiter=1,max_evaluations=2,
        verbose=false)

    # Regression: checking a STORED strict-portfolio path must evaluate
    # exactly the same household map as a fresh strict GE residual call.
    # Change only the first date's relative price slightly so that
    # foreign-bond revaluation makes the gross portfolio relevant.
    levctx=Stage3Context(short,copy(ctx.mu_pre),ctx.terms_pre,
                         1/16,0.0,true,"resume strict household portfolio")
    levinitial=copy(shortx)
    levinitial[2*Ttiny]=shortctx.terms_pre*1.001
    expected=similar(residual)
    actual=similar(residual)
    legacy=similar(residual)
    direct_work=TransitionWorkspace(short)
    resumed_work=TransitionWorkspace(short)
    legacy_work=TransitionWorkspace(short)
    FMR2020.stage3_flexible_full_residual!(
        expected,levinitial,levctx,direct_work;
        polish_households=true,strict_portfolio=true)
    economic_recheck_residual!(
        actual,levinitial,levctx,resumed_work;
        polish_households=true,strict_portfolio=true)
    economic_recheck_residual!(
        legacy,levinitial,levctx,legacy_work;
        polish_households=true,strict_portfolio=false)
    @test actual≈expected atol=1e-10 rtol=1e-10
    @test maximum(abs.(legacy.-expected))>1e-10
    @test all(isfinite,actual)
end
