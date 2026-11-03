
using Test
using FMR2020

@testset "Stage 3 shock and balance-sheet mechanisms" begin
    @testset "Unexpected reversal source timing" begin
        B=collect(1.0:1.0:200.0)
        result=credit_reversal_path(B,-4.0;shock_start=41,periods=200)
        @test length(result)==200
        @test result[1]==B[40]
        @test result[end]==-4.0
        @test isapprox(result[2],-4+(40+4)*(1-1/199)^4;atol=1e-12)
        @test_throws ArgumentError credit_reversal_path(B,-4.0;shock_start=1)
        # Refining the horizon must not redefine the original shock.
        oldshock=credit_reversal_path(B,-4.0;shock_start=41,periods=200)
        extension=credit_reversal_path(B,-4.0;shock_start=41,
                                       periods=300,decay_periods=200)
        @test extension[1:200]==oldshock
        @test all(extension[201:end] .== -4.0)
    end
    @testset "Positive trade-price roots with same-sign far endpoints" begin
        # This is the actual coefficient vector implied by the failed
        # T=200 initial aggregate guess in the owner's Windows log.
        C=-2.7782230139198565
        E=3.8257097540368417
        D=1.0149463596999195
        theta_star=3.0
        f(x)=C+E/x-D/x^theta_star
        xpeak=sqrt(3*D/E)
        @test f(0.001)<0
        @test f(1000.0)<0
        @test f(xpeak)>0
        root=FMR2020.foreign_debt_trade_root(C,E,D,theta_star)
        @test isapprox(root,0.7865895839646455;atol=1e-9)
        @test root<xpeak
        @test abs(f(root))<1e-9
        roots=FMR2020.foreign_debt_trade_roots(C,E,D,theta_star)
        @test isapprox(roots.low,root;atol=1e-11)
        @test roots.high !== nothing
        @test isapprox(roots.high,1.0379166453624395;atol=1e-9)
        @test roots.low<roots.peak<roots.high
        @test abs(f(roots.high))<1e-9
        # Opposite signs at two arbitrary far endpoints are NOT
        # required, but genuine lack of any positive root still fails.
        @test_throws DomainError FMR2020.foreign_debt_trade_root(
            -10.0,1.0,1.0,3.0)
    end

    @testset "Actual author MATLAB 399-vector has no first-period goods root" begin
        # Source: user-supplied exp_flex.mat and exp_flex_lev16.mat.
        # The author solved asset/labor residuals (<1.4e-8), but its
        # nested scalar fsolve(aggr,0.8) DID NOT solve t=1 goods clearing.
        # Do not silently call a non-root a valid nested trade-price solve.
        D=1.0149463620475911
        for (C,E,x_author,gap_author) in (
            (-2.7018950240416792,3.6621965903557725,
             0.9117778675311863,0.024336037042358205),
            (-2.696913415635473,3.6779032664825757,
             0.9098585114110708,0.002110419278378961))
            xstar=sqrt(3*D/E)
            goods_gap(x)=-C-E/x+D/x^3
            @test isapprox(goods_gap(x_author),gap_author;atol=1e-10)
            @test goods_gap(xstar)>0
            @test goods_gap(xstar)<=goods_gap(x_author)+1e-10
            @test_throws DomainError FMR2020.foreign_debt_trade_root(
                C,E,D,3.0)
            @test_throws DomainError FMR2020.foreign_debt_trade_roots(
                C,E,D,3.0)
        end
    end

    @testset "Regularized nonlinear direction remains descent near singular Jacobian" begin
        # The user's full GE stopped because the 399-variable Broyden
        # Jacobian was not refreshed after many near-zero damped steps.
        # The LM rescue direction must lower the local least-squares
        # objective even for a nearly singular finite-difference J.
        J=[1.0 1.0; 1.0 1.0+1e-11]
        residual=[0.02,-0.01]
        for λ in (1e-6,1e-4,1e-2,1.0,100.0)
            step=FMR2020.stage3_lm_direction(J,residual,λ)
            @test all(isfinite,step)
            grad=transpose(J)*residual
            @test sum(grad.*step)<0
        end
        @test_throws ArgumentError FMR2020.stage3_lm_direction(
            J,residual,0.0)
        @test_throws DimensionMismatch FMR2020.stage3_lm_direction(
            ones(2,3),residual,1e-3)
    end

    @testset "Unconverged checkpoint is never accepted as GE" begin
        mktempdir() do dir
            path=joinpath(dir,"stage3_unverified.chk")
            x=[1.0,2.0,3.0]
            res=[0.1,-0.03,0.02]
            FMR2020.stage3_write_checkpoint(path,x,res,2)
            @test isfile(path)
            @test FMR2020.stage3_read_checkpoint(path,3)==x
            @test_throws DimensionMismatch FMR2020.stage3_read_checkpoint(path,2)
            @test_throws DimensionMismatch FMR2020.stage3_write_checkpoint(
                path,x,[0.1],1)
            @test_throws ArgumentError FMR2020.stage3_write_checkpoint(
                path,x,res,0)
            @test_throws ArgumentError FMR2020.stage3_write_checkpoint(
                path,[NaN,2.0],[0.0,0.0],1)
        end
    end

    @testset "MATLAB full-grid and constrained EGM reference" begin
        params=default_parameters()
        gr=make_asset_grids(params,baseline_reference(params).K)
        layout=source_matlab_grid_layout(gr)
        @test layout.coarse_total==193
        @test layout.fine_total==769
        @test layout.coarse_borrow_index==33
        @test layout.fine_borrow_index==129
        @test layout.coarse_below==32
        @test layout.fine_below==128
        @test layout.coarse_active==161
        @test layout.fine_active==641
        @test length(gr.coarse)==161
        @test length(gr.coarse_full)==193
        @test gr.coarse_full[33]==gr.borrowing_limit
        chain=rouwenhorst(params)
        full_future=repeat(reshape(1.0 .+ 0.02 .* max.(gr.coarse_full,0.0),
                                   :,1),1,params.n_s)
        source=source_matlab_egm_step(
            gr,chain,params,full_future,0.98,0.01,0.01,
            1.0,1.8,1.0,0.25)
        @test size(source.coarse_consumption)==(193,7)
        @test size(source.fine_consumption)==(769,7)
        # MATLAB interp1 without extrapolation can return NaN for
        # off-grid nodes. In its vectorized masked sum, 0*NaN is NaN.
        # Do NOT misreport these as valid household decisions.
        valid=isfinite.(source.fine_consumption)
        @test any(valid)
        @test all(source.fine_consumption[valid] .> 0)
        @test all(isnan.(source.fine_savings) .== .!valid)
        @test all(source.fine_savings[valid] .>= gr.borrowing_limit)
        # Enforced lower saving bound in the constrained region follows
        # author backsolve_egm.m, not a blanket household-debt clip.
        @test all(isfinite,source.endogenous_threshold)
        # The source constraint binds iff the current asset a is
        # BELOW the state-specific endogenous cutoff; no common
        # binding status at a=bmin is implied across productivity.
        for j in 1:params.n_s
            binding=findall(a->a<=source.endogenous_threshold[j],
                            gr.fine_full)
            occupied_valid=[i for i in binding
                            if isfinite(source.fine_savings[i,j])]
            @test all(source.fine_savings[occupied_valid,j] .==
                      gr.borrowing_limit)
        end
        @test_throws DimensionMismatch source_matlab_egm_step(
            gr,chain,params,full_future[33:end,:],0.98,0.01,0.01,
            1.0,1.8,1.0,0.25)
    end

    @testset "Zero exchange-rate shock preserves household net wealth" begin
        grid=Float64[-2,-1,0,1,2,3,4]
        mu=zeros(7,2)
        mu[3,1]=0.1
        mu[4,1]=0.3
        mu[5,2]=0.6
        dest=similar(mu)
        info=revalue_foreign_debt!(dest,mu,grid,2.0,1.0,-0.5,
                                    1.2,1.2,1/16)
        @test abs(info.mass_after-1)<1e-12
        @test abs(info.net_wealth_before-info.net_wealth_after)<1e-12
        @test maximum(abs.(dest.-mu))<1e-12
        @test isapprox(info.capital+info.bonds,1.5;atol=1e-12)
        @test_throws ArgumentError revalue_foreign_debt!(
            mu,mu,grid,2.0,1.0,-0.5,1.2,1.2,1/16)
    end
    @testset "Strict aggregate gross capital and foreign bond conservation" begin
        grid=Float64[-3,-2,-1,0,1,2,3,4]
        mu=zeros(8,1)
        mu[3,1]=0.2 # pre-shock net wealth -1
        mu[5,1]=0.2 # net wealth +1
        mu[7,1]=0.6 # net wealth +3
        # Mass-weighted pre-shock net wealth: 1.8 = q0*K0 + Bprev.
        legacy=revalue_foreign_debt!(
            similar(mu),mu,grid,2.0,1.0,-0.2,1.0,1.2,1/16)
        dest=similar(mu)
        strict=revalue_foreign_debt!(
            dest,mu,grid,2.0,1.0,-0.2,1.0,1.2,1/16;
            strict_portfolio=true)
        @test abs(legacy.capital_market_gap)>1e-3
        @test isapprox(strict.capital,2.0;atol=1e-12)
        @test isapprox(strict.bonds,-0.2;atol=1e-12)
        @test abs(strict.capital_market_gap)<1e-12
        @test abs(strict.inherited_bond_market_gap)<1e-12
        @test isapprox(strict.net_wealth_before,1.8;atol=1e-12)
        @test isapprox(strict.net_wealth_after,2.0-0.2/1.2;
                       atol=1e-12)
        @test isapprox(sum(dest),1.0;atol=1e-12)
        @test strict.capital_normalization>1.0
        @test !legacy.strict_portfolio
        # Strict conservation correctly rejects inconsistent inherited
        # total net wealth, rather than silently manufacturing capital.
        @test_throws DomainError revalue_foreign_debt!(
            similar(mu),mu,grid,2.0,1.0,0.0,1.0,1.2,1/16;
            strict_portfolio=true)
    end
    @testset "Exchange-rate change revalues liabilities and conserves mass" begin
        grid=Float64[-3,-2,-1,0,1,2,3,4]
        mu=zeros(8,1)
        mu[3,1]=0.5  # borrower with negative foreign bond holdings
        mu[6,1]=0.5
        dest=similar(mu)
        info=revalue_foreign_debt!(dest,mu,grid,2.0,1.0,0.0,
                                    1.0,2.0,0.0)
        @test isapprox(sum(dest),1.0;atol=1e-12)
        @test isapprox(info.net_wealth_after,
                       info.capital+info.bonds/2;atol=1e-12)
        @test minimum(dest)>=0
        @test isapprox(info.net_wealth_after-info.net_wealth_before,
                       0.25;atol=1e-12)
    end
    @testset "Positive probability leaving full asset grid errors" begin
        grid=Float64[-1,0,1]
        mu=zeros(3,1)
        mu[1,1]=1.0
        dest=similar(mu)
        @test_throws DomainError revalue_foreign_debt!(
            dest,mu,grid,1.0,1.0,0.0,1.0,0.1,0.0)
    end
end

@testset "Explicit 599-variable flexible source-independent closure" begin
    # 199 K plus 200 L plus 200 terms: no scalar nested fsolve.
    params=default_parameters()
    steady=solve_three_steady_states(params;verbose=false)
    ctx=stage3_expansion_context(steady)
    w=FMR2020.TransitionWorkspace(ctx.economy)
    n=length(ctx.economy.B)
    @test n==200
    extended=FMR2020.stage3_extend_horizon(ctx,240)
    @test extended.economy.B[1:200]==ctx.economy.B
    @test all(extended.economy.B[201:end] .== ctx.economy.terminal.B)
    @test isapprox(ctx.economy.B[end],ctx.economy.terminal.B;
                   atol=1e-10,rtol=0)
    @test extended.mu_pre==ctx.mu_pre
    x0=vcat(FMR2020.initial_transition_guess(ctx.economy),ones(n))
    longer=FMR2020.stage3_prolong_aggregate_guess(x0,ctx,extended)
    @test length(longer)==3*240-1
    @test longer[1:199]==x0[1:199]
    @test longer[240:439]==x0[200:399]
    @test longer[480:679]==x0[400:599]
    @test_throws ArgumentError FMR2020.stage3_prolong_aggregate_guess(
        x0,ctx,extended;rho=1.0)
    @test length(x0)==599
    @test_throws DimensionMismatch FMR2020.stage3_flexible_full_prices!(
        w,ctx,x0[1:end-1])
    FMR2020.stage3_flexible_full_prices!(w,ctx,x0)
    markup=(params.epsilon_f-1)/params.epsilon_f
    for t in 1:n
        @test isapprox(w.wage[t]*w.L[t]/
            ((1-params.alpha)*w.Y[t]),markup;atol=1e-12)
        @test isapprox(w.rental[t]*w.K[t]/
            (params.alpha*w.Y[t]),markup;atol=1e-12)
    end
    # The first two blocks are the unchanged original Stage3 equations.
    fill!(w.aggregate_savings,0.0)
    fill!(w.aggregate_home,1.0)
    full=zeros(3*n-1)
    core=zeros(2*n-1)
    FMR2020.stage3_market_residual!(core,w,ctx)
    FMR2020.stage3_flexible_full_market_residual!(full,w,ctx)
    @test full[1:2*n-1]==core
    # Independently impose the paper's nonlinear union wage condition.
    # The source MATLAB incorrectly uses L=RHS when varphi != 1.
    p=ctx.economy.params
    ratio=p.chi*(p.epsilon_w-1)/p.epsilon_w*
          w.wage[1]/(ctx.economy.reference.chi_dis*w.aggregate_home[1])
    @test isapprox(w.labor_implied[1],ratio^(1/p.varphi);atol=1e-12)
    @test all(isfinite,full)
    # T is the stationary continuation boundary, not a final date
    # where the source artificially copies consumption from date T-1.
    # Its 2 residuals pin L[T] and terms[T] to the calibrated SS.
    base=copy(full)
    w.aggregate_home[n]=1.25
    FMR2020.stage3_flexible_full_market_residual!(full,w,ctx)
    @test full==base
    old_L=w.L[n]
    w.L[n]+=0.25
    FMR2020.stage3_flexible_full_market_residual!(full,w,ctx)
    scaleL=max(1.0,ctx.economy.initial.L)
    @test isapprox(full[2*n-1]-base[2*n-1],0.25/scaleL;atol=1e-12)
    w.L[n]=old_L
    old_terms=w.terms[n]
    w.terms[n]+=0.1
    FMR2020.stage3_flexible_full_market_residual!(full,w,ctx)
    @test isapprox(full[end]-base[end],
         0.1/max(1.0,abs(ctx.economy.terminal.terms));atol=1e-12)
    fixed_a=zeros(3*n-1)
    fixed_b=zeros(3*n-1)
    FMR2020.stage3_market_residual!(fixed_a,w,ctx;fixed=true)
    w.terms[n]=old_terms
    FMR2020.stage3_market_residual!(fixed_b,w,ctx;fixed=true)
    @test isapprox(fixed_a[end]-fixed_b[end],
         0.1/max(1.0,abs(ctx.economy.terminal.terms));atol=1e-12)
    @test_throws DimensionMismatch FMR2020.stage3_flexible_full_market_residual!(
        zeros(399),w,ctx)
end


@testset "Economic fixed-FX solver: independent acceptance boundaries" begin
    # Include the real command-line solver without executing a full GE root.
    # This also verifies imports and methods used by the command-line script.
    include(joinpath(@__DIR__,"..","scripts","solve_economic_fixed_ge.jl"))
    good=(maximum_scaled_solver_residual=7e-7,
          max_asset_market_level=2e-6,
          max_labor_market_level=1e-6,
          max_goods_market_level=3e-6,
          max_household_budget=1e-14,
          max_mass_deviation=1e-14,
          max_nonlinear_phillips=1e-12,
          max_terminal_household_budget=1e-8,
          terminal_K_gap=1e-7,
          terminal_L_gap=1e-7,
          terminal_terms_gap=1e-8,
          max_dynamic_euler=5e-5,
          weighted_mean_dynamic_euler=1e-6)
    @test fixed_market_acceptance(good)
    @test fixed_euler_acceptance(good)
    @test fixed_terminal_acceptance(good,20.0,1.0)
    @test !fixed_terminal_acceptance(
        merge(good,(max_terminal_household_budget=0.001,)),20.0,1.0)
    @test !fixed_terminal_acceptance(
        merge(good,(terminal_K_gap=0.005,)),20.0,1.0)
    @test fixed_market_acceptance(merge(good,(max_dynamic_euler=3.8e-4,)))
    @test !fixed_euler_acceptance(merge(good,(max_dynamic_euler=3.8e-4,)))
    @test !fixed_market_acceptance(merge(good,(max_goods_market_level=1e-2,)))
    @test !fixed_market_acceptance(merge(good,(max_nonlinear_phillips=1e-3,)))
    @test !fixed_market_acceptance(merge(good,(max_mass_deviation=1e-5,)))
    mktempdir() do d
        path=joinpath(d,"shorter.csv")
        fields="period,B,K,L,Y,I,q,r,w,terms,foreign_consumption,home_consumption,wealth,nominal_exchange_rate"
        open(path,"w") do io
            println(io,fields)
            for t in 1:40
                println(io,join((t,0.0,1.5,1.0,2.0,0.1,1.0,0.02,
                                0.5,1.0,0.1,1.0,1.0,1.0),","))
            end
        end
        longer=(K=fill(1.5,40),L=fill(1.0,40),
                Y=fill(2.0,40),I=fill(0.1,40),
                wage=fill(0.5,40),rate=fill(0.02,40),
                terms=fill(1.0,40),aggregate_home=fill(1.0,40))
        comparison=fixed_horizon_comparison(path,longer)
        @test comparison.maximum==0.0
        @test comparison.first_periods==40
        longer.K[5]+=0.03
        disturbed=fixed_horizon_comparison(path,longer)
        @test isapprox(disturbed.by_variable["K"],0.03/1.5;atol=1e-12)
        @test isapprox(disturbed.maximum,0.03/1.5;atol=1e-12)
    end
end


@testset "Independent one-state Euler/KKT best-response root" begin
    # Deterministic logarithmic-consumption benchmark with affine future
    # policy c_next(a')=1+a'; beta*R_next=1 and current resources=2.
    # Exact interior solution: 1/(2-a') = 1/(1+a'), hence a'=0.5.
    grid=Float64[0.0,1.0,2.0]
    future=reshape(Float64[1.0,2.0,3.0],3,1)
    prob=Float64[1.0]
    best=FMR2020.household_best_response_saving(
        grid,future,prob,1.0,1.0,1.0,2.0,1.0,0.0)
    @test !best.binding
    @test isapprox(best.savings,0.5;atol=1e-9)
    @test isapprox(best.consumption,1.5;atol=1e-9)
    @test abs(best.euler_surplus)<1e-9

    # The unconstrained root would be negative; borrowing KKT must bind.
    constrained=FMR2020.household_best_response_saving(
        grid,future,prob,1.0,1.0,1.0,0.5,1.0,0.0)
    @test constrained.binding
    @test constrained.savings==0.0
    @test constrained.consumption==0.5
    @test constrained.euler_surplus>=0.0

    # Never silently clamp optimal saving at the numerical upper grid bound.
    @test_throws DomainError FMR2020.household_best_response_saving(
        grid,future,prob,1.0,1.0,1.0,7.0,1.0,0.0)
    @test_throws DimensionMismatch FMR2020.household_best_response_saving(
        grid,future,Float64[0.4,0.6],1.0,1.0,1.0,2.0,1.0,0.0)
    @test_throws ArgumentError FMR2020.household_best_response_saving(
        grid,future,Float64[0.99],1.0,1.0,1.0,2.0,1.0,0.0)

    # Coarse EGM interpolation can miss the continuous Euler optimum;
    # the slope-reusing local corrector must recover the *same* independent
    # 0.5 savings root without altering the household budget.
    slopes=similar(future)
    FMR2020.pchip_slopes!(view(slopes,:,1),grid,view(future,:,1))
    corrected=FMR2020.euler_refine_one_saving(
        grid,future,slopes,prob,1.0,1.0,1.0,2.0,1.0,0.0,0.48)
    @test corrected.refined
    @test !corrected.offgrid
    @test isapprox(corrected.saving,0.5;atol=1e-8)
    @test isapprox(corrected.consumption,1.5;atol=1e-8)
    @test abs(corrected.relative_gap)<1e-9
    @test isapprox(corrected.saving+corrected.consumption,2.0;atol=1e-12)
    unchanged=FMR2020.euler_refine_one_saving(
        grid,future,slopes,prob,1.0,1.0,1.0,2.0,1.0,0.0,0.5)
    @test !unchanged.refined
    @test !unchanged.offgrid
    bound=FMR2020.euler_refine_one_saving(
        grid,future,slopes,prob,1.0,1.0,1.0,0.5,1.0,0.0,0.1)
    @test bound.refined
    @test bound.saving==0.0
    @test isapprox(bound.consumption,0.5;atol=1e-12)
    # Off-grid constrained roots are not coerced into a fake solution.
    outside=FMR2020.euler_refine_one_saving(
        grid,future,slopes,prob,1.0,1.0,1.0,7.0,1.0,0.0,1.0)
    @test outside.offgrid
    @test outside.saving==1.0
end
