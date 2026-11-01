#!/usr/bin/env julia
# Nirei (2015) Issue #5 Stage 8: diagnostic-only finite-N investment wedge.
# Same frozen firm state and same iid future-productivity vector in EVERY
# component. Conditioning is explicitly on existence of positive-X fixed
# point. Source model, ES1, and original iid draw stream are unchanged.
#
# julia --project=. scripts/diagnostic_investment_wedge.jl [T] [N] [seed] [draws] [max_attempts]
# Default T=501 includes frozen dates 0,4,20,100,500.
include(joinpath(@__DIR__, "diagnostic_shock_rejection.jl"))
using Random: randn!
using Statistics: median, quantile, cov

const WEDGE_DATES = Set((0, 4, 20, 100, 500))
const WEDGE_DRAW_TAG = 0x243f6a8885a308d3

# The production workspace already contains all 2N thresholds and the
# power-weight a_next^(rho/alpha) * d^rho / N after enumeration.
# K_fixed is Appendix A's exp-transformed expected log capital, NOT
# necessarily the level implied by the continuum Eq. (18). The latter
# remains the reference for the original production ES1 selection.
function fixed_capital_response(state::SimulationState, K_fixed::Float64)
    ws = state.workspace
    firms = state.firms
    pars = state.calibration.parameters
    rho = state.calibration.constants.rho
    n = length(firms.capital)
    x = 0.0
    s = 0.0
    @inbounds for i in 1:n
        lambda_i = firms.lambda[i]
        kappa = K_fixed <= ws.switches[2i-1] ? -1 :
                K_fixed <= ws.switches[2i] ? 0 : 1
        d = ws.depreciated[i]
        multiplier = lambda_i^kappa
        x += d*(multiplier-1.0)/n
        s += ws.base_power[i]*lambda_i^(rho*kappa)
    end
    return (X=x, Phi=s^(inv(rho)), power=s)
end

# Evaluate the exact signed investment change from K_fixed to K_ES1.
# These are threshold-crossing identities, not a derivative of a step
# function and not a separate counterfactual general equilibrium.
function crossing_accounting(state::SimulationState,
                             K_fixed::Float64, K_ES1::Float64)
    ws=state.workspace
    firms=state.firms
    rho=state.calibration.constants.rho
    n=length(firms.capital)
    dx=0.0
    dpower=0.0
    nchanged=0
    nlow=0
    nhigh=0
    @inbounds for i in 1:n
        lo=ws.switches[2i-1]
        hi=ws.switches[2i]
        old_kappa=K_fixed<=lo ? -1 : K_fixed<=hi ? 0 : 1
        new_kappa=K_ES1<=lo ? -1 : K_ES1<=hi ? 0 : 1
        if old_kappa!=new_kappa
            nchanged+=1
            lambda_i=firms.lambda[i]
            dx+=ws.depreciated[i]*
                (lambda_i^new_kappa-lambda_i^old_kappa)/n
            dpower+=ws.base_power[i]*
                (lambda_i^(rho*new_kappa)-lambda_i^(rho*old_kappa))
            a=min(old_kappa,new_kappa)
            b=max(old_kappa,new_kappa)
            nlow+=(a==-1 && b>=0)
            nhigh+=(a<=0 && b==1)
        end
    end
    return (delta_X=dx, delta_power=dpower,
            changed_firms=nchanged, lower_events=nlow,
            upper_events=nhigh, total_events=nlow+nhigh)
end

function _wedge_mean_se(values::Vector{Float64})
    n=length(values)
    n>0 || return Dict("n"=>0,"mean"=>"undefined","SE"=>"undefined")
    if n==1
        return Dict("n"=>1,"mean"=>only(values),"SE"=>"undefined")
    end
    return Dict(
        "n"=>n,"mean"=>mean(values),"SE"=>std(values)/sqrt(n),
        "median"=>median(values),"sd"=>std(values),
        "p05"=>quantile(values,0.05),
        "p95"=>quantile(values,0.95),
        "p_negative"=>count(x->x<0.0,values)/n
    )
end

# The same independent MC seed construction used by the previous 32-draw
# capital-closure diagnostic; no change to the production/retry RNGs.
function wedge_frozen_state!(state::SimulationState, date::Int,
                             draws::Int, seed::Int)
    state.period==date || error("frozen date does not match state.period")
    pars=state.calibration.parameters
    c=state.calibration.constants
    re=state.re
    prior_K=state.K
    prior_mu=state.mu_prior
    prior_capital=copy(state.firms.capital)
    prior_productivity=copy(state.firms.productivity)
    !state.firms.next_ready || error("unexpected pending shock on frozen state")
    expected=compute_expected_equilibrium(re,state.calibration,prior_K,prior_mu)
    reference_K18=(1.0-pars.delta)*prior_K+
        state.calibration.threshold.A^(c.rho-1.0)/c.rho*expected.X_next
    expected_K=expected.K_next
    base=UInt64(seed) ⊻ (UInt64(date+1)*0x9e3779b97f4a7c15) ⊻ WEDGE_DRAW_TAG
    rng=Xoshiro(base)
    # Verify independent draw reproducibility before ANY fixed-point work.
    test_rng=Xoshiro(base)
    check=similar(state.firms.next_productivity)
    randn!(test_rng,check)
    @inbounds for i in eachindex(check)
        check[i]=exp(c.mean_log_a+pars.sigma_log_a*check[i])
    end
    rows=Vector{Dict{String,Any}}()
    feasible=0
    fixed_X_all=Float64[]
    fixed_X18_all=Float64[]
    x_fixed=Float64[]
    x_fixed18=Float64[]
    x_es1=Float64[]
    wedge_forecast=Float64[]
    wedge_equilibrium=Float64[]
    wedge_total=Float64[]
    wedge_forecast18=Float64[]
    wedge_equilibrium18=Float64[]
    logeps=Float64[]
    ranks=Int[]
    max_selection_gaps=Float64[]
    max_identity_residual=0.0
    max_crossing_residual=0.0
    max_CES_residual=0.0
    max_reference_crossing_residual=0.0
    max_reference_CES_residual=0.0
    for draw in 1:draws
        next_a=state.firms.next_productivity
        randn!(rng,next_a)
        @inbounds @simd for i in eachindex(next_a)
            next_a[i]=exp(c.mean_log_a+pars.sigma_log_a*next_a[i])
        end
        if draw==1
            next_a==check || error("independent MC seed/draw is not deterministic")
        end
        state.firms.next_ready=true
        fixed_points=enumerate_fixed_points!(state.workspace,state.firms,
                                             state.calibration)
        response=fixed_capital_response(state,expected_K)
        response18=fixed_capital_response(state,reference_K18)
        reference_cross=crossing_accounting(state,expected_K,reference_K18)
        max_reference_crossing_residual=max(
            max_reference_crossing_residual,
            abs(reference_cross.delta_X-(response18.X-response.X))
        )
        max_reference_CES_residual=max(
            max_reference_CES_residual,
            abs(reference_cross.delta_power-(response18.power-response.power))
        )
        isapprox(reference_cross.delta_X,response18.X-response.X;
                 atol=1e-9,rtol=1e-9) ||
            error("Eq18-vs-AppendixA fixed-K threshold X identity failed")
        isapprox(reference_cross.delta_power,response18.power-response.power;
                 atol=1e-9,rtol=1e-9) ||
            error("Eq18-vs-AppendixA fixed-K CES identity failed")
        push!(fixed_X_all,response.X)
        push!(fixed_X18_all,response18.X)
        ws=state.workspace
        nfix=fixed_points.total
        npositive=fixed_points.admissible
        row=Dict{String,Any}(
            "t"=>date,"draw"=>draw,"feasible"=>npositive>0,
            "X_expected"=>expected.X_next,"K_expected_AppendixA"=>expected_K,
            "K_expected_Eq18"=>reference_K18,
            "K_expected_discrepancy"=>reference_K18-expected_K,
            "X_fixed"=>response.X,"Phi_fixed"=>response.Phi,
            "X_fixed_Eq18"=>response18.X,
            "Phi_fixed_Eq18"=>response18.Phi,
            "reference_fixed_X_difference"=>response18.X-response.X,
            "reference_changed_firms"=>reference_cross.changed_firms,
            "reaction_residual_at_K_expected"=>response.Phi-expected_K,
            "n_fixed"=>nfix,"n_positive_X"=>npositive,
            "min_fixed_K"=>fixed_points.minimum_K,
            "max_fixed_K"=>fixed_points.maximum_K,
            "closest_fixed_abs_logK_distance"=>
                minimum(abs(log(ws.fixed_capital[j]/expected_K)) for j in 1:nfix)
        )
        if npositive>0
            feasible+=1
            choice=select_es1(ws,expected.X_next,reference_K18)
            j=choice.index
            xe=ws.fixed_investment[j]
            ke=ws.fixed_capital[j]
            cross=crossing_accounting(state,expected_K,ke)
            fw=response.X-expected.X_next
            ew=xe-response.X
            tw=xe-expected.X_next
            fw18=response18.X-expected.X_next
            ew18=xe-response18.X
            abs(fw18+ew18-tw)<1e-11 ||
                error("Eq18 reference wedge decomposition failed")
            power_change=ke^c.rho-response.power
            cross_res=abs(cross.delta_X-ew)
            ces_res=abs(cross.delta_power-power_change)
            wedge_res=abs(fw+ew-tw)
            max_identity_residual=max(max_identity_residual,wedge_res)
            max_crossing_residual=max(max_crossing_residual,cross_res)
            max_CES_residual=max(max_CES_residual,ces_res)
            isapprox(cross.delta_X,ew;atol=1e-9,rtol=1e-9) ||
                error("threshold crossings fail to explain equilibrium X wedge")
            isapprox(cross.delta_power,power_change;atol=1e-9,rtol=1e-9) ||
                error("threshold crossings fail to explain CES power difference")
            abs(wedge_res)<1e-11 || error("investment wedge identity failed")
            pos_idxs=[i for i in 1:nfix if ws.fixed_investment[i]>0.0]
            max_pos_x=maximum(ws.fixed_investment[i] for i in pos_idxs)
            min_pos_x=minimum(ws.fixed_investment[i] for i in pos_idxs)
            closest_pos_K=minimum(abs(log(ws.fixed_capital[i]/expected_K))
                                   for i in pos_idxs)
            row["X_ES1"]=xe
            row["K_ES1"]=ke
            row["forecast_wedge"]=fw
            row["equilibrium_wedge"]=ew
            row["total_wedge"]=tw
            row["forecast_wedge_Eq18"]=fw18
            row["equilibrium_wedge_Eq18"]=ew18
            row["epsilon"]=log(xe/expected.X_next)
            row["ES1_capital_rank"]=j
            row["ES1_rank_fraction"]=j/nfix
            row["ES1_minus_max_positive_X"]=xe-max_pos_x
            row["min_positive_X"]=min_pos_x
            row["max_positive_X"]=max_pos_x
            row["closest_positive_abs_logK_distance"]=closest_pos_K
            row["ES1_abs_logK_distance"]=abs(log(ke/expected_K))
            row["changed_firms"]=cross.changed_firms
            row["lower_switch_events"]=cross.lower_events
            row["upper_switch_events"]=cross.upper_events
            row["threshold_events"]=cross.total_events
            row["delta_X_from_crossings"]=cross.delta_X
            row["delta_CES_power_from_crossings"]=cross.delta_power
            push!(x_fixed,response.X)
            push!(x_fixed18,response18.X)
            push!(x_es1,xe)
            push!(wedge_forecast,fw)
            push!(wedge_equilibrium,ew)
            push!(wedge_total,tw)
            push!(wedge_forecast18,fw18)
            push!(wedge_equilibrium18,ew18)
            push!(logeps,row["epsilon"])
            push!(ranks,j)
            push!(max_selection_gaps,xe-max_pos_x)
        end
        push!(rows,row)
    end
    state.firms.next_ready=false
    state.period==date && state.K==prior_K && state.mu_prior==prior_mu &&
    state.firms.capital==prior_capital &&
    state.firms.productivity==prior_productivity ||
        error("frozen MC mutated the economy's predetermined state")
    # Exact conditional means use the SAME feasible draws for every term.
    summary=Dict{String,Any}(
        "t"=>date,"draws"=>draws,"n_feasible"=>feasible,
        "n_infeasible"=>draws-feasible,
        "feasibility_probability_MC"=>feasible/draws,
        "K_expected_AppendixA"=>expected_K,
        "K_reference_Eq18_for_ES1"=>reference_K18,
        "X_expected_AppendixA"=>expected.X_next,
        "X_fixed_unconditional_iid"=>_wedge_mean_se(fixed_X_all),
        "X_fixed_Eq18_unconditional_iid"=>_wedge_mean_se(fixed_X18_all),
        "X_fixed_conditional_feasible"=>_wedge_mean_se(x_fixed),
        "X_fixed_Eq18_conditional_feasible"=>_wedge_mean_se(x_fixed18),
        "X_ES1_conditional_feasible"=>_wedge_mean_se(x_es1),
        "forecast_wedge_conditional"=>_wedge_mean_se(wedge_forecast),
        "equilibrium_wedge_conditional"=>_wedge_mean_se(wedge_equilibrium),
        "total_wedge_conditional"=>_wedge_mean_se(wedge_total),
        "forecast_wedge_conditional_Eq18"=>_wedge_mean_se(wedge_forecast18),
        "equilibrium_wedge_conditional_Eq18"=>_wedge_mean_se(wedge_equilibrium18),
        "paired_equilibrium_minus_forecast_wedge"=>_wedge_mean_se(
            wedge_equilibrium .- wedge_forecast
        ),
        "epsilon_conditional_feasible"=>_wedge_mean_se(logeps),
        "mean_ES1_rank"=>isempty(ranks) ? "undefined" : mean(ranks),
        "share_ES1_is_max_positive_X"=>
            isempty(max_selection_gaps) ? "undefined" :
            count(x->abs(x)<1e-11,max_selection_gaps)/length(max_selection_gaps),
        "max_abs_wedge_identity_residual"=>max_identity_residual,
        "max_abs_threshold_X_residual"=>max_crossing_residual,
        "max_abs_threshold_CES_power_residual"=>max_CES_residual,
        "max_abs_reference_X_residual"=>max_reference_crossing_residual,
        "max_abs_reference_CES_residual"=>max_reference_CES_residual,
        "covariance_forecast_equilibrium_wedges"=>
            length(wedge_forecast)>1 ?
                cov(wedge_forecast,wedge_equilibrium) : "undefined",
        "conditioning"=>"All three wedge means condition on identical iid draws with at least one positive-X fixed point; nonexistence is not zero investment.",
        "fixed_K_response_is_not_equilibrium"=>true,
        "original_iid_unconditional_epsilon_undefined_on_failures"=>true
    )
    if feasible>0
        f=mean(wedge_forecast)
        e=mean(wedge_equilibrium)
        total=mean(wedge_total)
        abs(f+e-total)<1e-10 ||
            error("conditional wedge means do not sum")
        summary["mean_wedge_additivity_residual"]=f+e-total
        f18=mean(wedge_forecast18)
        e18=mean(wedge_equilibrium18)
        abs(f18+e18-total)<1e-10 ||
            error("Eq18 reference conditional wedge means do not sum")
        summary["mean_wedge_Eq18_additivity_residual"]=f18+e18-total
        if abs(total)>1e-12
            summary["forecast_Eq18_signed_share"]=f18/total
            summary["equilibrium_Eq18_signed_share"]=e18/total
            summary["forecast_signed_share"]=f/total
            summary["equilibrium_signed_share"]=e/total
        end
    end
    return (summary=summary,rows=rows)
end

function write_wedge_csv(path::String, rows::Vector{Dict{String,Any}})
    cols=[
        "t","draw","feasible","X_expected","K_expected_AppendixA",
        "K_expected_Eq18","K_expected_discrepancy","X_fixed","Phi_fixed",
        "X_fixed_Eq18","Phi_fixed_Eq18","reference_fixed_X_difference",
        "reference_changed_firms",
        "reaction_residual_at_K_expected","n_fixed","n_positive_X",
        "min_fixed_K","max_fixed_K","closest_fixed_abs_logK_distance",
        "X_ES1","K_ES1","forecast_wedge","equilibrium_wedge",
        "total_wedge","forecast_wedge_Eq18",
        "equilibrium_wedge_Eq18","epsilon","ES1_capital_rank","ES1_rank_fraction",
        "ES1_minus_max_positive_X","min_positive_X","max_positive_X",
        "closest_positive_abs_logK_distance","ES1_abs_logK_distance",
        "changed_firms","lower_switch_events","upper_switch_events",
        "threshold_events","delta_X_from_crossings",
        "delta_CES_power_from_crossings"
    ]
    open(path,"w") do io
        println(io,join(cols,","))
        for row in rows
            println(io,join((string(get(row,key,"")) for key in cols),","))
        end
    end
end

function run_stage8_wedge()
    T=length(ARGS)>=1 ? parse(Int,ARGS[1]) : 501
    N=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 110000
    seed=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 20261101
    draws=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 32
    max_attempts=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 100
    T>=5 && N>0 && draws>=2 && max_attempts>=1 ||
        error("invalid arguments")
    state=initialize_simulation(;N=N,base_seed=seed)
    retry_stream=Xoshiro(UInt64(seed) ⊻ RETRY_SEED_TAG)
    allrows=Vector{Dict{String,Any}}()
    summaries=Vector{Dict{String,Any}}()
    rejected=0
    for t in 0:(T-1)
        if t in WEDGE_DATES
            sample=wedge_frozen_state!(state,t,draws,seed)
            append!(allrows,sample.rows)
            push!(summaries,sample.summary)
            println("WEDGE DATE=",t," summary=",sample.summary)
        end
        out=one_diagnostic_period!(state,retry_stream;max_attempts=max_attempts)
        out.failed && error("state reconstruction failure at t=$t: "*out.failure)
        rejected+=out.rejected
        if N==110000 && seed==20261101 && t==0
            isapprox(out.record.K_next_firm,8.546731346350889;
                     atol=1e-9,rtol=1e-11) ||
                error("original first-period ES1 regression failure")
            isapprox(out.record.X,0.06936827009543808;
                     atol=1e-10,rtol=1e-11) ||
                error("original first-period investment regression failure")
        end
    end
    root=normpath(joinpath(@__DIR__,"..","results"))
    mkpath(root)
    stem="diagnostic_investment_wedge_T$(T)_N$(N)_seed$(seed)"
    write_wedge_csv(joinpath(root,stem*".csv"),allrows)
    summary=Dict{String,Any}(
        "experiment"=>"issue5_stage8_frozen_state_investment_wedge",
        "not_published_replication"=>true,
        "configuration"=>Dict(
            "T"=>T,"N"=>N,"seed"=>seed,"draws_per_state"=>draws,
            "max_attempts"=>max_attempts,"julia_version"=>string(VERSION)
        ),
        "status"=>Dict(
            "completed_periods"=>T,"total_rejections_on_reconstructed_path"=>rejected,
            "n_frozen_states"=>length(summaries)
        ),
        "conditional_MC"=>summaries,
        "interpretation"=>"Frozen-state iid Monte Carlo. ES1 outcomes conditional on positive-investment equilibrium existence. Path reconstruction after t=5 uses diagnostic shock rejection."
    )
    open(joinpath(root,stem*".toml"),"w") do io
        TOML.print(io,summary)
    end
    println("WEDGE COMPLETED states=",length(summaries),
            " rows=",length(allrows)," reconstructed_dates=",T,
            " rejections=",rejected)
    println("WEDGE FILES ",stem*".csv | ",stem*".toml")
    return nothing
end

if abspath(PROGRAM_FILE)==@__FILE__
    run_stage8_wedge()
end
