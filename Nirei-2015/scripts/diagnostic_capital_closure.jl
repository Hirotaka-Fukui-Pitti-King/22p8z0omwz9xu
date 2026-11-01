#!/usr/bin/env julia
# Issue #5: finite-N CES capital accounting and rational-expectations
# closure audit along the already defined shock-rejection diagnostic path.
# PRODUCTION src/ and ES1 UNCHANGED. Shock rejection changes shock law.
#
# julia --project=. scripts/diagnostic_capital_closure.jl [T] [N] [seed] [max_attempts]
#
# This script includes the existing diagnostic rejection simulation so it
# uses the EXACT SAME draw/ES1/update code (no separate firm transition).
include(joinpath(@__DIR__, "diagnostic_shock_rejection.jl"))
using Random: randn!
using Statistics: median

const DETAIL_DATES=Set((0,1,2,3,4,5,10,20,50,100,250,500,750,1000,1099))
const MC_DATES=Set((0,4,20,100,500))
const MC_DRAWS=32

function diagnostic_capital_identity(
    state, initial_firm_K::Vector{Float64}, period_record::PeriodRecord
)
    cap=state.calibration
    pars=cap.parameters
    rho=cap.constants.rho
    alpha=pars.alpha
    countfirms=length(initial_firm_K)
    countfirms==length(state.firms.capital) ||
        error("Firm count changed over this date")
    capital_after=state.firms.capital
    productivity_after=state.firms.productivity
    # Both choices of the CES power sum are evaluated directly from the
    # FIRM arrays rather than from previously computed aggregate residuals.
    base_S=0.0
    policy_plus=0.0
    policy_minus=0.0
    total_S=0.0
    total_X=0.0
    nplus=0
    nminus=0
    nzero=0
    @inbounds for i in eachindex(initial_firm_K)
        d=(1.0-pars.delta)*initial_firm_K[i]
        ka=capital_after[i]
        lambda_i=state.firms.lambda[i]
        a=productivity_after[i]
        isfinite(d) && d>0 && isfinite(a) && a>0 && ka>0 ||
            error("Capital/productivity not positive at firm $i")
        unadjusted=a^(rho/alpha)*d^rho
        realized=a^(rho/alpha)*ka^rho
        base_S+=unadjusted
        total_S+=realized
        total_X+=ka-d
        ratio=ka/d
        if isapprox(ratio,lambda_i;rtol=1e-11,atol=1e-13)
            nplus+=1
            policy_plus+=realized-unadjusted
        elseif isapprox(ratio,inv(lambda_i);rtol=1e-11,atol=1e-13)
            nminus+=1
            policy_minus+=realized-unadjusted
        elseif isapprox(ratio,1.0;rtol=1e-11,atol=1e-13)
            nzero+=1
            abs(realized-unadjusted)<1e-10*max(1.0,abs(realized)) ||
                error("Inactive firm has nonzero CES adjustment")
        else
            error("Updated firm capital is not one of {-1,0,1} choices at $i")
        end
    end
    S0=base_S/countfirms
    Splus=policy_plus/countfirms
    Sminus=policy_minus/countfirms
    S1=total_S/countfirms
    Sp=Splus+Sminus
    check_power=S0+Sp-S1
    eq_power=period_record.K_next_firm^rho
    mean_X=total_X/countfirms
    isapprox(S0+Sp,S1;rtol=2e-11,atol=1e-12) ||
        error("Exact CES power decomposition failed")
    isapprox(S1,eq_power;rtol=2e-11,atol=1e-12) ||
        error("CES power sum differs from firm's actual GE capital")
    isapprox(mean_X,period_record.X;rtol=2e-11,atol=1e-12) ||
        error("Firm average investment differs from aggregate investment")
    nplus==period_record.n_investing &&
    nminus==period_record.n_divesting &&
    nzero==period_record.n_inactive ||
        error("Firm {-1,0,+1} choice counts do not match production record")
    oldK=period_record.K_prior
    Knew=period_record.K_next_firm
    growth=log(Knew/oldK)
    noadjust_log=log(S0/oldK^rho)/rho
    policy_log=log1p(Sp/S0)/rho
    resid=growth-noadjust_log-policy_log
    abs(resid)<1e-9 ||
        error("Exact log capital accounting failed: "*string(resid))
    factor=cap.threshold.A^(rho-1.0)/rho
    K18=(1.0-pars.delta)*oldK+factor*period_record.X
    K18>0.0 || error("Continuum equation (18) predicts nonpositive K")
    pred_loggrowth=log(K18/oldK)
    # The linear Eq (26) is already present in PeriodRecord and can
    # also be independently reconstructed from its recorded forecast.
    linear_pred=period_record.expected_K_next
    logK_expect=log(linear_pred/cap.steady_state.K)
    xinnov=log(period_record.X/period_record.expected_X)
    linear_logdev=logK_expect+pars.delta*xinnov
    linear_res=log(Knew/cap.steady_state.K)-linear_logdev
    isapprox(linear_res,period_record.capital_linear_gap;
             rtol=1e-10,atol=1e-10) ||
        error("Independent Equation (26) residual differs from production")
    logmu=log(period_record.mu/cap.steady_state.mu)
    tautological_cont_res=logmu-
        state.re.h*log(Knew/cap.steady_state.K)
    abs(tautological_cont_res)<1e-10 ||
        error("Diagnostic realized mu deviates from implemented closure")
    return Dict{String,Any}(
        "t"=>period_record.t,
        "K_prior"=>oldK,
        "K_next"=>Knew,
        "log_K_growth"=>growth,
        "X_actual"=>period_record.X,
        "X_expected_Appendix_A"=>period_record.expected_X,
        "epsilon"=>period_record.epsilon,
        "S0_no_investment"=>S0,
        "Splus_investment"=>Splus,
        "Sminus_divestment"=>Sminus,
        "Sp_net_policy"=>Sp,
        "S_total_direct"=>S1,
        "max_CES_power_identity_residual"=>max(
            abs(check_power),abs(S1-eq_power)
        ),
        "exact_log_no_adjust_contribution"=>noadjust_log,
        "exact_log_policy_contribution"=>policy_log,
        "exact_log_accounting_residual"=>resid,
        "equation18_implied_K"=>K18,
        "equation18_capital_residual"=>Knew-K18,
        "equation18_loggrowth_residual"=>growth-pred_loggrowth,
        "eq26_linear_log_capital_residual"=>linear_res,
        "mu_realized"=>period_record.mu,
        "mu_expected"=>period_record.expected_mu,
        "logK_expectation_error"=>period_record.logK_expectation_error,
        "logmu_expectation_error"=>period_record.logmu_expectation_error,
        "realized_mu_stable_path_tautological_residual"=>
            tautological_cont_res,
        "n_investing"=>nplus,
        "n_divesting"=>nminus,
        "n_inactive"=>nzero,
        "n_fixed"=>period_record.n_fixed,
        "n_positive_X"=>period_record.n_positive_X
    )
end

function _independent_conditional_mc!(
    state, t::Int; draws::Int=MC_DRAWS, seed::Int=20261101
)
    state.period==t || error("Conditional MC state date mismatch")
    cap=state.calibration
    p=cap.parameters
    c=cap.constants
    forecast=compute_expected_equilibrium(
        state.re,cap,state.K,state.mu_prior
    )
    expected_K=(1-p.delta)*state.K+
        cap.threshold.A^(c.rho-1)/c.rho*forecast.X_next
    rng=Xoshiro(UInt64(seed) ⊻
        (UInt64(t+1)*0x9e3779b97f4a7c15) ⊻ 0x243f6a8885a308d3)
    preK=state.K
    premu=state.mu_prior
    precap=copy(state.firms.capital)
    prea=copy(state.firms.productivity)
    periods_before=state.period
    xvalid=Float64[]
    evalid=Float64[]
    logKvalid=Float64[]
    logmuvalid=Float64[]
    bad=0
    for _ in 1:draws
        nexta=state.firms.next_productivity
        randn!(rng,nexta)
        @inbounds @simd for i in eachindex(nexta)
            nexta[i]=exp(c.mean_log_a+p.sigma_log_a*nexta[i])
        end
        state.firms.next_ready=true
        result=enumerate_fixed_points!(state.workspace,state.firms,cap)
        if result.admissible==0
            bad+=1
        else
            j=select_es1(state.workspace,forecast.X_next,expected_K).index
            x=state.workspace.fixed_investment[j]
            k=state.workspace.fixed_capital[j]
            push!(xvalid,x)
            push!(evalid,log(x/forecast.X_next))
            push!(logKvalid,log(k/cap.steady_state.K))
            # This is the CONDITIONAL mean of the CURRENT Julia
            # realized-state closure, not a proven RE shadow price.
            push!(logmuvalid,state.re.h*log(k/cap.steady_state.K))
        end
    end
    state.firms.next_ready=false
    state.K==preK && state.mu_prior==premu &&
    state.period==periods_before &&
    state.firms.capital==precap &&
    state.firms.productivity==prea ||
        error("Conditional Monte Carlo modified predetermined firm/macro state")
    function avg_se(values::Vector{Float64})
        if isempty(values)
            return ("undefined","undefined")
        elseif length(values)==1
            return (values[1],"undefined: n=1")
        end
        return (mean(values),std(values)/sqrt(length(values)))
    end
    xmean,xse=avg_se(xvalid)
    emean,ese=avg_se(evalid)
    kmean,kse=avg_se(logKvalid)
    mumean,muse=avg_se(logmuvalid)
    return Dict{String,Any}(
        "date"=>t,
        "MC_draws"=>draws,
        "no_positive_X_count"=>bad,
        "positive_X_equilibrium_count"=>length(xvalid),
        "nonexistence_frequency"=>bad/draws,
        "mean_X_conditional_on_positive_X"=>xmean,
        "SE_X_conditional"=>xse,
        "mean_epsilon_conditional_on_positive_X"=>emean,
        "SE_epsilon_conditional"=>ese,
        "mean_logK_conditional_on_positive_X"=>kmean,
        "SE_logK_conditional"=>kse,
        "mean_logmu_under_implemented_closure_conditional"=>mumean,
        "SE_logmu_conditional"=>muse,
        "appendixA_expected_logK"=>forecast.logK_next,
        "appendixA_expected_logmu"=>forecast.logmu_next,
        "appendixA_expected_X"=>forecast.X_next,
        "conditional_logK_minus_AppendixA_forecast"=>
            isempty(logKvalid) ? "undefined" : kmean-forecast.logK_next,
        "conditional_logmu_minus_AppendixA_forecast"=>
            isempty(logmuvalid) ? "undefined" : mumean-forecast.logmu_next,
        "conditioning"=>"positive-X equilibrium exists; draws without X>0 excluded from conditional means, not counted as zero",
        "closure_not_validated_by_conditional_MC"=>true
    )
end

function _format_csv(path::String, rows::Vector{Dict{String,Any}})
    isempty(rows) && return
    cols=sort(collect(keys(first(rows))))
    open(path,"w") do io
        println(io,join(cols,","))
        for row in rows
            println(io,join((string(get(row,k,"")) for k in cols),","))
        end
    end
end

function run_capital_closure_audit()
    horizon=length(ARGS)>=1 ? parse(Int,ARGS[1]) : 1100
    n=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 110000
    seed=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 20261101
    max_attempts=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 100
    horizon>=6 && n>=1 && max_attempts>=1 ||
        throw(ArgumentError("invalid simulation parameters"))
    state=initialize_simulation(;N=n,base_seed=seed)
    original_benchmark=initialize_simulation(;N=n,base_seed=seed)
    retry_stream=Xoshiro(UInt64(seed) ⊻ RETRY_SEED_TAG)
    calc=state.calibration
    re=state.re
    # Sensitivities calculated DIRECTLY from the paper's B^{-1}D.
    # log X^e_{t+1} = (F11 log K_{t+1} + F12 log mu_t
    #                   - (1-delta)log K_{t+1})/delta.
    dlogXe_dlogmu=re.F[1,2]/calc.parameters.delta
    dlogXe_dlogK_stable=(
        re.F[1,1]+re.F[1,2]*re.h-
        (1-calc.parameters.delta)
    )/calc.parameters.delta
    rows=Vector{Dict{String,Any}}()
    recs=PeriodRecord[]
    conditional=Vector{Dict{String,Any}}()
    rejected_total=0
    rejected_first100=0
    reject_dates=0
    failure="none"
    for t in 0:(horizon-1)
        # Conditional future-shock MC is independent from ORIGINAL shocks
        # and sampled before the date's realized shock is generated.
        if t in MC_DATES
            push!(conditional,_independent_conditional_mc!(
                state,t;draws=MC_DRAWS,seed=seed
            ))
        end
        cap_before=copy(state.firms.capital)
        out=one_diagnostic_period!(
            state,retry_stream;max_attempts=max_attempts
        )
        rejected_total+=out.rejected
        rejected_first100+=t<100 ? out.rejected : 0
        reject_dates+=out.rejected>0 ? 1 : 0
        if out.failed
            failure=out.failure
            println("CAPITAL AUDIT FAILED at t=",t,": ",failure)
            break
        end
        rec=out.record::PeriodRecord
        if n==110000 && seed==20261101 && t<=4
            verify_initial_prefix!(state,original_benchmark,rec,t)
        end
        row=diagnostic_capital_identity(state,cap_before,rec)
        row["rejected_shocks"]=out.rejected
        push!(rows,row)
        push!(recs,rec)
        if t in DETAIL_DATES || (t+1)%100==0
            println("CAPITAL AUDIT t=",t,
                " K=",rec.K_next_firm,
                " gK=",row["log_K_growth"],
                " no_adjust=",row["exact_log_no_adjust_contribution"],
                " policy=",row["exact_log_policy_contribution"],
                " eq18_res=",row["equation18_capital_residual"],
                " lin_res=",row["eq26_linear_log_capital_residual"],
                " epsilon=",rec.epsilon,
                " rejects=",out.rejected)
        end
    end
    root=normpath(joinpath(@__DIR__,"..","results"))
    mkpath(root)
    stem="diagnostic_capital_closure_T"*string(horizon)*
         "_N"*string(n)*"_seed"*string(seed)
    full=joinpath(root,stem*".csv")
    points=joinpath(root,stem*"_key_dates.csv")
    toml=joinpath(root,stem*".toml")
    _format_csv(full,rows)
    _format_csv(points,filter(r->r["t"] in DETAIL_DATES,rows))
    totalg=sum((r["log_K_growth"] for r in rows);init=0.0)
    noadj=sum((r["exact_log_no_adjust_contribution"] for r in rows);init=0.0)
    policy=sum((r["exact_log_policy_contribution"] for r in rows);init=0.0)
    init_K=isempty(rows) ? state.K : first(rows)["K_prior"]
    actualcumulative=isempty(rows) ? 0.0 :
        log(last(rows)["K_next"]/init_K)
    abs(totalg-noadj-policy)<1e-7 ||
        error("Cumulative exact CES decomposition failed")
    abs(totalg-actualcumulative)<1e-7 ||
        error("Cumulative log growth not equal to endpoint log ratio")
    post=filter(r->r["t"]>=100,rows)
    mean_over(which,key)=isempty(which) ? 0.0 :
        mean(r[key] for r in which)
    max_over(which,key)=isempty(which) ? 0.0 :
        maximum(abs(r[key]) for r in which)
    structural_source_audit=Dict{String,Any}(
        "equation18_derivation"=>"Appendix B integrates investments over STATIONARY UNIFORM firm gaps and bounded log productivity increments; it is not literally an identity for finite-N realized CES capital with arbitrary gap distributions or Gaussian log productivity.",
        "paper_expectations"=>"Proposition 5 proves local saddle determinacy for Appendix A FIRST-ORDER EXPECTATION system, not a global finite-N stochastic equilibrium with positive-X existence everywhere.",
        "replication_mu_closure"=>"src/simulation.jl imposes realized logmu_t = h * logK_(t+1); by construction the realized stable-jump residual is zero and cannot itself validate that assumption.",
        "no_paper_justified_replacement_implemented"=>true,
        "shock_rejection_changed_conditional_shock_law"=>true
    )
    summary=Dict{String,Any}(
        "experiment"=>"diagnostic_capital_collapse_closure",
        "is_published_replication"=>false,
        "configuration"=>Dict(
            "N"=>n,"seed"=>seed,"periods"=>horizon,
            "max_attempts"=>max_attempts,"julia_version"=>string(VERSION),
            "conditional_MC_dates"=>sort(collect(MC_DATES)),
            "conditional_MC_draws_each"=>MC_DRAWS
        ),
        "status"=>Dict(
            "completed_periods"=>length(rows),
            "completed"=>length(rows)==horizon && failure=="none",
            "failure"=>failure,
            "rejected_profiles"=>rejected_total,
            "rejected_first100"=>rejected_first100,
            "rejected_after100"=>rejected_total-rejected_first100,
            "dates_with_rejection"=>reject_dates
        ),
        "exact_decomposition"=>Dict(
            "cumulative_log_growth"=>totalg,
            "cumulative_log_noadjust_contribution"=>noadj,
            "cumulative_log_policy_contribution"=>policy,
            "endpoint_log_growth"=>actualcumulative,
            "endpoint_K"=>isempty(rows) ? init_K :
                last(rows)["K_next"],
            "initial_K"=>init_K,
            "max_per_date_accounting_residual"=>
                max_over(rows,"exact_log_accounting_residual"),
            "max_CES_power_sum_residual"=>
                max_over(rows,"max_CES_power_identity_residual"),
            "mean_gK"=>mean_over(rows,"log_K_growth"),
            "mean_noadjust_log_contribution"=>
                mean_over(rows,"exact_log_no_adjust_contribution"),
            "mean_policy_log_contribution"=>
                mean_over(rows,"exact_log_policy_contribution")
        ),
        "approximation_residuals"=>Dict(
            "mean_eq18_capital_residual"=>mean_over(rows,
                "equation18_capital_residual"),
            "max_abs_eq18_capital_residual"=>max_over(rows,
                "equation18_capital_residual"),
            "mean_linear_Eq26_log_residual"=>mean_over(rows,
                "eq26_linear_log_capital_residual"),
            "max_abs_linear_Eq26_log_residual"=>max_over(rows,
                "eq26_linear_log_capital_residual"),
            "mean_eps"=>mean_over(rows,"epsilon"),
            "median_eps"=>isempty(rows) ? 0.0 :
                median(r["epsilon"] for r in rows)
        ),
        "post100"=>Dict(
            "mean_gK"=>mean_over(post,"log_K_growth"),
            "mean_eps"=>mean_over(post,"epsilon"),
            "mean_eq18_capital_residual"=>mean_over(post,
                "equation18_capital_residual"),
            "mean_linear_Eq26_log_residual"=>mean_over(post,
                "eq26_linear_log_capital_residual")
        ),
        "macro_feedback_coefficients"=>Dict(
            "h_realized_continuation"=>re.h,
            "partial_dlogExpectedX_dlogPriorMu_at_fixedK"=>
                dlogXe_dlogmu,
            "dlogExpectedX_dlogPriorK_along_realized_stable_path"=>
                dlogXe_dlogK_stable,
            "F11"=>re.F[1,1],"F12"=>re.F[1,2]
        ),
        "conditional_acceptance_MC"=>conditional,
        "source_audit"=>structural_source_audit,
        "outputs"=>Dict(
            "all_periods_csv"=>full,
            "selected_dates_csv"=>points
        )
    )
    open(toml,"w") do io
        TOML.print(io,summary)
    end
    println("CAPITAL AUDIT RESULTS ",summary["exact_decomposition"])
    println("CAPITAL AUDIT LINEAR ",summary["approximation_residuals"])
    println("CAPITAL AUDIT FEEDBACK ",summary["macro_feedback_coefficients"])
    println("CAPITAL AUDIT CONDITIONAL MC ",conditional)
    println("CAPITAL AUDIT SUMMARY completed=",length(rows),"/",horizon,
        " rejects=",rejected_total," status=",failure)
    println("CAPITAL AUDIT FILES ",full," | ",points," | ",toml)
    length(rows)==horizon && failure=="none" ||
        error("Capital closure audit stopped after "*string(length(rows))*
              " periods: "*failure)
    return nothing
end

if abspath(PROGRAM_FILE)==@__FILE__
    run_capital_closure_audit()
end
