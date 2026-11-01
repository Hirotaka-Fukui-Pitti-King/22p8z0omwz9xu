#!/usr/bin/env julia
# Diagnostic-only finite-N Nirei (2015) shock rejection experiment, Issue #5.
# From Nirei-2015:
# julia --project=. scripts/diagnostic_shock_rejection.jl [T] [N] [seed] [max_attempts]
# This is NOT the published Appendix C simulation: it samples productivity
# conditional on the EXISTENCE of a positive-X equilibrium. Never overwrite
# src/simulation.jl, production ES1, or the original benchmark output paths.

using Nirei2015
using Random: Xoshiro
using Statistics: mean, std, cor
using TOML

const EXPERIMENT_LABEL = "diagnostic_shock_rejection_counterfactual"
const RETRY_SEED_TAG = 0x6a09e667f3bcc909

function trial_profile(
    state,candidate_expected_K::Float64; accepted::Bool=false
)
    firms=state.firms
    ws=state.workspace
    n=length(firms.capital)
    cap=state.calibration
    rho=cap.constants.rho
    exponent=rho/cap.parameters.alpha
    sum_loga=0.0
    sum_loga2=0.0
    sum_CES=0.0
    investing=0
    divesting=0
    @inbounds for i in 1:n
        a=accepted ? firms.productivity[i] : firms.next_productivity[i]
        la=log(a)
        sum_loga+=la
        sum_loga2+=la^2
        sum_CES+=a^exponent*ws.depreciated[i]^rho
        investing+=(candidate_expected_K>ws.switches[2i])
        divesting+=(candidate_expected_K<=ws.switches[2i-1])
    end
    mean_log=sum_loga/n
    return (
        mean_log_a=mean_log,
        variance_log_a=max(0.0,sum_loga2/n-mean_log^2),
        CES_baseline_power=sum_CES/n,
        invest_share_at_expected_K=investing/n,
        divest_share_at_expected_K=divesting/n
    )
end

# Rejected realizations must not advance any predetermined state.
function check_unmodified_state(
    state, oldK::Float64, oldmu::Float64, oldt::Int,
    before_capital::Vector{Float64}, before_productivity::Vector{Float64}
)
    state.K==oldK && state.mu_prior==oldmu && state.period==oldt ||
        error("Rejected shock modified macro state or simulation date")
    state.firms.capital==before_capital ||
        error("Rejected shock modified firm capital")
    state.firms.productivity==before_productivity ||
        error("Rejected shock modified current firm productivity")
    state.firms.next_ready ||
        error("Rejected shock did not leave a pending productivity profile")
    return nothing
end

# Retry RNG is separated from the production stream and does not consume
# the next date's FIRST-attempt benchmark shock. FirmRNGs is immutable,
# so we temporarily substitute a new struct containing the SAME original
# first three RNGs and an independent subsequent-productivity RNG.
function retry_rngs(original::FirmRNGs, rng::Xoshiro)
    return FirmRNGs(
        original.lambda_seed,
        original.initial_productivity_seed,
        original.initial_gap_seed,
        original.subsequent_productivity_seed,
        original.lambda_rng,
        original.initial_productivity_rng,
        original.initial_gap_rng,
        rng
    )
end

function is_only_positive_equilibrium_failure(state, failure::SimulationFailure)
    ws=state.workspace
    n=ws.n_fixed
    return n>0 && state.firms.next_ready &&
        occursin("No equilibrium has X>0", failure.reason) &&
        all(x->x<=0.0, view(ws.fixed_investment,1:n))
end

function one_diagnostic_period!(
    state, retry_stream::Xoshiro; max_attempts::Int=100
)
    max_attempts>=1 || throw(ArgumentError("max_attempts must be positive"))
    original_rngs=state.rngs
    t=state.period
    oldK=state.K
    oldmu=state.mu_prior
    before_capital=copy(state.firms.capital)
    before_productivity=copy(state.firms.productivity)
    expected=compute_expected_equilibrium(state.re,state.calibration,oldK,oldmu)
    p=state.calibration.parameters
    rho=state.calibration.constants.rho
    expected_K=(1-p.delta)*oldK+
        state.calibration.threshold.A^(rho-1)/rho*expected.X_next
    attempts=Vector{Dict{String,Any}}()
    result=nothing
    failure="none"
    for attempt in 1:max_attempts
        # ONLY attempt 1 consumes the original benchmark stream.
        state.rngs= attempt==1 ? original_rngs :
            retry_rngs(original_rngs,retry_stream)
        try
            rec=simulate_period!(state)
            # simulate_period! calls the ORIGINAL exact enumerator, ES1,
            # firm transition, CES check, and realized macro closure.
            result=rec
            accepted_profile=trial_profile(state,expected_K;accepted=true)
            # The latest accepted productivity vector was moved into
            # firms.productivity by production's advance_productivity!.
            push!(attempts,Dict{String,Any}(
                "attempt"=>attempt,"accepted"=>true,
                "n_fixed"=>rec.n_fixed,"n_positive_X"=>rec.n_positive_X,
                "expected_X"=>rec.expected_X,
                "actual_X"=>rec.X,
                "epsilon"=>rec.epsilon,
                "mean_log_a"=>accepted_profile.mean_log_a,
                "variance_log_a"=>accepted_profile.variance_log_a,
                "CES_baseline_power"=>accepted_profile.CES_baseline_power,
                "invest_share_at_expected_K"=>
                    accepted_profile.invest_share_at_expected_K,
                "divest_share_at_expected_K"=>
                    accepted_profile.divest_share_at_expected_K
            ))
            break
        catch caught
            if caught isa SimulationFailure &&
               is_only_positive_equilibrium_failure(state,caught)
                check_unmodified_state(
                    state,oldK,oldmu,t,before_capital,before_productivity
                )
                ws=state.workspace
                profile=trial_profile(state,expected_K)
                push!(attempts,Dict{String,Any}(
                    "attempt"=>attempt,"accepted"=>false,
                    "n_fixed"=>ws.n_fixed,
                    "n_positive_X"=>0,
                    "max_equilibrium_X"=>
                        maximum(view(ws.fixed_investment,1:ws.n_fixed)),
                    "mean_log_a"=>profile.mean_log_a,
                    "variance_log_a"=>profile.variance_log_a,
                    "CES_baseline_power"=>profile.CES_baseline_power,
                    "invest_share_at_expected_K"=>
                        profile.invest_share_at_expected_K,
                    "divest_share_at_expected_K"=>
                        profile.divest_share_at_expected_K
                ))
            else
                failure=sprint(showerror,caught)
                break
            end
        finally
            # All future first attempts use the benchmark stream.
            state.rngs=original_rngs
        end
    end
    if result===nothing && failure=="none"
        failure="No positive-X equilibrium within maximum attempted profiles"
    end
    return (
        record=result,
        draws=attempts,
        rejected=count(x->!x["accepted"],attempts),
        failed=failure!="none",
        failure=failure,
        prior_K=oldK,
        prior_mu=oldmu
    )
end

function lagone(v::Vector{Float64})
    return length(v)>=3 && std(v)>0 ?
        cor(view(v,1:length(v)-1),view(v,2:length(v))) : 0.0
end

function moments_table2(records::Vector{PeriodRecord}, burnin::Int)
    length(records)>burnin+2 || return Dict{String,Any}()
    retained=records[(burnin+1):end]
    keys=("logY","logC","logX","logL","logK")
    values=Dict{String,Vector{Float64}}(
        k=>[Float64(getproperty(r,Symbol(k))) for r in retained]
        for k in keys
    )
    correlations=Dict{String,Float64}()
    for i in 1:length(keys), j in (i+1):length(keys)
        correlations[keys[i]*"__"*keys[j]]=cor(values[keys[i]],values[keys[j]])
    end
    return Dict{String,Any}(
        "retained"=>length(retained),
        "sd_log_levels"=>Dict(k=>std(values[k]) for k in keys),
        "lag1_autocorrelation"=>Dict(k=>lagone(values[k]) for k in keys),
        "cross_correlations"=>correlations,
        "mean_X"=>mean(r.X for r in retained),
        "mean_expected_X"=>mean(r.expected_X for r in retained),
        "mean_epsilon"=>mean(r.epsilon for r in retained),
        "sd_epsilon"=>std([r.epsilon for r in retained]),
        "lag1_epsilon"=>lagone([r.epsilon for r in retained])
    )
end

function emit_csv(path::String, records::Vector{PeriodRecord})
    open(path,"w") do io
        fields=fieldnames(PeriodRecord)
        println(io,join(string.(fields),","))
        for rec in records
            println(io,join((string(getfield(rec,k)) for k in fields),","))
        end
    end
end

function emit_attempt_csv(path::String, rows::Vector{Dict{String,Any}})
    names=(
        "t","attempt","accepted","n_fixed","n_positive_X",
        "max_equilibrium_X","mean_log_a","variance_log_a",
        "CES_baseline_power","invest_share_at_expected_K",
        "divest_share_at_expected_K","expected_X","actual_X","epsilon"
    )
    open(path,"w") do io
        println(io,join(names,","))
        for row in rows
            println(io,join((string(get(row,k,"")) for k in names),","))
        end
    end
end

function verify_initial_prefix!(state, benchmark, record::PeriodRecord, t::Int)
    if t<=4
        expected=simulate_period!(benchmark)
        record.K_next_firm==expected.K_next_firm &&
        record.X==expected.X &&
        record.expected_X==expected.expected_X ||
            error("Diagnostic first-attempt deviates from benchmark at t=$t")
        state.firms.capital==benchmark.firms.capital &&
        state.firms.productivity==benchmark.firms.productivity ||
            error("Firm state differs from benchmark in pre-rejection prefix")
    end
end

function selection_effect_summaries(
    attempts::Vector{Dict{String,Any}}
)
    vars=(
        "mean_log_a","variance_log_a","CES_baseline_power",
        "invest_share_at_expected_K","divest_share_at_expected_K"
    )
    accepted=filter(x->x["accepted"],attempts)
    rejected=filter(x->!x["accepted"],attempts)
    observed=Dict{String,Any}(
        "accepted_draws"=>length(accepted),
        "rejected_draws"=>length(rejected)
    )
    for key in vars
        observed[key*"_accepted_mean"]=
            isempty(accepted) ? "not defined" :
            mean(x[key] for x in accepted)
        observed[key*"_rejected_mean"]=
            isempty(rejected) ? "not defined" :
            mean(x[key] for x in rejected)
    end
    return observed
end

function run_experiment()
    T=length(ARGS)>=1 ? parse(Int,ARGS[1]) : 20
    N=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 110000
    seed=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 20261101
    max_attempts=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 100
    T>=2 && N>0 && seed>=0 && max_attempts>=1 ||
        throw(ArgumentError("invalid diagnostic configuration"))
    state=initialize_simulation(;N=N,base_seed=seed)
    benchmark=initialize_simulation(;N=N,base_seed=seed)
    retry_stream=Xoshiro(UInt64(seed) ⊻ RETRY_SEED_TAG)
    initial_K=state.K
    records=PeriodRecord[]
    attempt_rows=Vector{Dict{String,Any}}()
    period_rows=Vector{Dict{String,Any}}()
    prefix_check= N==110000 && seed==20261101
    failure="none"
    for t in 0:(T-1)
        outcome=one_diagnostic_period!(
            state,retry_stream;max_attempts=max_attempts
        )
        for a in outcome.draws
            row=copy(a)
            row["t"]=t
            push!(attempt_rows,row)
        end
        rejected=outcome.rejected
        push!(period_rows,Dict{String,Any}(
            "t"=>t,
            "attempt_count"=>length(outcome.draws),
            "rejected_count"=>rejected,
            "prior_K"=>outcome.prior_K,
            "prior_mu"=>outcome.prior_mu,
            "status"=>outcome.failed ? "failure" : "accepted"
        ))
        if outcome.failed
            failure=outcome.failure
            println("DIAGNOSTIC FAILURE t=",t,
                "; rejections=",rejected,"; cause=",failure)
            break
        end
        rec=outcome.record::PeriodRecord
        if prefix_check
            verify_initial_prefix!(state,benchmark,rec,t)
            if t==5
                # The primary t=5 shock must have failed exactly as in the
                # original benchmark; the accepted retry then advances.
                rejected>=1 ||
                    error("Diagnostic did not reject known t=5 original shock")
            end
        end
        push!(records,rec)
        if t<=9 || (t+1)%100==0 || rejected>0
            println("DIAGNOSTIC t=",t,
                " attempts=",length(outcome.draws),
                " rejected=",rejected,
                " X=",rec.X," expected_X=",rec.expected_X,
                " epsilon=",rec.epsilon," K=",rec.K_next_firm)
        end
    end
    root=normpath(joinpath(@__DIR__,"..","results"))
    mkpath(root)
    stem=EXPERIMENT_LABEL*"_T"*string(T)*"_N"*string(N)*"_seed"*string(seed)
    csv=joinpath(root,stem*".csv")
    attempts_csv=joinpath(root,stem*"_attempts.csv")
    toml=joinpath(root,stem*".toml")
    emit_csv(csv,records)
    emit_attempt_csv(attempts_csv,attempt_rows)
    nrej=count(x->!x["accepted"],attempt_rows)
    accepted=length(records)
    summary=Dict{String,Any}(
        "experiment"=>EXPERIMENT_LABEL,
        "is_published_replication"=>false,
        "original_iid_shocks_not_preserved_after_rejection"=>true,
        "julia_version"=>string(VERSION),
        "configuration"=>Dict(
            "firms"=>N,"periods"=>T,"seed"=>seed,
            "max_attempts"=>max_attempts,"burnin"=>min(100,T-1),
            "initial_K"=>initial_K,
            "retry_RNG"=>"independent fixed Xoshiro stream",
            "first_attempt_rng"=>"original Nirei2015 production stream"
        ),
        "status"=>Dict(
            "completed"=>accepted==T && failure=="none",
            "completed_periods"=>accepted,
            "requested_periods"=>T,
            "failure"=>failure
        ),
        "rejections"=>Dict(
            "total_rejected_profiles"=>nrej,
            "total_attempts"=>length(attempt_rows),
            "accepted_periods"=>accepted,
            "fraction_dates_with_rejection"=>
                isempty(period_rows) ? 0.0 :
                count(x->x["rejected_count"]>0,period_rows)/length(period_rows),
            "maximum_attempts_in_one_period"=>
                isempty(period_rows) ? 0 :
                maximum(x->x["attempt_count"],period_rows),
            "mean_attempts_per_completed_period"=>
                accepted==0 ? 0.0 : length(attempt_rows)/accepted,
            "per_period"=>period_rows
        ),
        "integrity"=>Dict(
            "max_abs_fixed_point_residual"=>
                isempty(records) ? 0.0 :
                maximum(x->abs(x.fixed_point_residual),records),
            "max_abs_investment_residual"=>
                isempty(records) ? 0.0 :
                maximum(x->abs(x.investment_residual),records),
            "max_abs_stable_jump_residual"=>
                isempty(records) ? 0.0 :
                maximum(x->abs(x.stable_jump_residual),records),
            "baseline_prefix_exact_check"=>
                prefix_check && accepted>=5,
        ),
        "post_burnin"=>moments_table2(records,min(100,T-1)),
        "accepted_vs_rejected_productivity"=>selection_effect_summaries(
            attempt_rows
        ),
        "rejection_blocks"=>Dict(
            "initial_100_rejections"=>
                sum(x["rejected_count"] for x in period_rows if x["t"]<100),
            "subsequent_rejections"=>
                sum(x["rejected_count"] for x in period_rows if x["t"]>=100)
        ),
        "output"=>Dict(
            "period_csv"=>csv,"attempt_csv"=>attempts_csv,
            "summary_toml"=>toml
        ),
        "interpretation"=>"Diagnostic rejection samples productivity conditional on positive-X equilibrium existence, not author's verified original code or a self-consistent new RE equilibrium."
    )
    open(toml,"w") do io
        TOML.print(io,summary)
    end
    println("DIAGNOSTIC SUMMARY completed=",accepted,"/",T,
        " rejected=",nrej," total_draws=",length(attempt_rows),
        " status=",failure)
    println("DIAGNOSTIC FILES ",toml," ; ",csv," ; ",attempts_csv)
    if !isempty(summary["post_burnin"])
        post=summary["post_burnin"]
        println("DIAGNOSTIC POST-BURNIN ",post)
        println("DIAGNOSTIC SELECTION EFFECTS ",
            summary["accepted_vs_rejected_productivity"])
        println("DIAGNOSTIC REJECTION BLOCKS ",summary["rejection_blocks"])
        println("DIAGNOSTIC INITIAL K=",initial_K,
            " FINAL K=",last(records).K_next_firm,
            " FINAL X=",last(records).X)
    end
    accepted==T && failure=="none" ||
        error("Diagnostic shock rejection stopped after $accepted periods: $failure")
    return nothing
end

if abspath(PROGRAM_FILE)==@__FILE__
    run_experiment()
end
