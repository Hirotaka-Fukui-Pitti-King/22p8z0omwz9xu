#!/usr/bin/env julia
# Issue #5 Stage 13: DIAGNOSTIC-only dynamic J=2 choice-set counterfactual.
# Original src/, calibration, ES1, and productivity RNG stream are untouched.
#
# julia --project=. scripts/diagnostic_multistep_simulation.jl smoke
# julia --project=. scripts/diagnostic_multistep_simulation.jl run 20 110000 20261101
# julia --project=. scripts/diagnostic_multistep_simulation.jl run 100 110000 20261101
# julia --project=. scripts/diagnostic_multistep_simulation.jl run 1100 110000 20261101
#
# run T executes BOTH A (J=1 for t<5, J=2 afterward) and B (J=2 always).
# Each experiment starts from an independent clone of the same RNG seed.
# A verified equilibrium/domain failure is reported in TOML, not bypassed.
include(joinpath(@__DIR__, "diagnostic_multistep_equilibrium.jl"))
using Statistics: mean, std

const STAGE13_REFERENCE_K5 = 8.145907725689481
const STAGE13_REFERENCE_X5 = 0.16603788812361434

function stage13_period_row(rec::PeriodRecord, n::Int)
    return Dict{String,Any}(
        "t"=>rec.t, "J"=>1,"K_prior"=>rec.K_prior,
        "K_next"=>rec.K_next_firm,"X"=>rec.X,"Y"=>rec.Y,"C"=>rec.C,
        "L"=>rec.L,"mu"=>rec.mu,"epsilon"=>rec.epsilon,
        "expected_K"=>rec.expected_K_next,"expected_X"=>rec.expected_X,
        "n_fixed"=>rec.n_fixed,"n_positive_X"=>rec.n_positive_X,
        "invest_1"=>rec.n_investing,"invest_2plus"=>0,
        "divest_1"=>rec.n_divesting,"divest_2plus"=>0,
        "inactive"=>rec.n_inactive,
        "CES_residual"=>rec.fixed_point_residual,
        "investment_residual"=>rec.investment_residual,
        "linear_K_gap"=>rec.capital_linear_gap,
        "logK_forecast_error"=>rec.logK_expectation_error,
        "logmu_forecast_error"=>rec.logmu_expectation_error,
        "stable_jump_residual"=>rec.stable_jump_residual,
        "level_resource_gap"=>(rec.Y-rec.C-rec.X)/rec.Y,
        "linear_resource_residual"=>rec.resources_linear_residual,
        "linear_euler_residual"=>rec.euler_linear_residual,
        "firm_count"=>n,
        "shock_checksum"=>"" # t=5 checksum is populated in multistep step.
    )
end

function stage13_multistep_step!(state::SimulationState, J::Int)
    J>=2 || error("Stage 13 multi-step step requires J>=2")
    t=state.period
    kp=state.K
    mup=state.mu_prior
    cal=state.calibration
    p=cal.parameters
    steady=cal.steady_state
    rho=cal.constants.rho
    expected=compute_expected_equilibrium(state.re,cal,kp,mup)
    draw_productivity!(state.firms,state.rngs,cal)
    a_next=state.firms.next_productivity
    shock_checksum=sum(log, a_next)
    enumeration=multistep_enumerate(state.firms,cal,J)
    positives=0
    best=0
    objective=Inf
    capital_tiebreak=Inf
    ke18=(1-p.delta)*kp+cal.threshold.A^(rho-1)/rho*expected.X_next
    for (idx,root) in enumerate(enumeration.roots)
        if root.X>0.0
            positives+=1
            score=abs(log(root.X/expected.X_next))
            distance=abs(log(root.K/ke18))
            if score<objective ||
                (score==objective && distance<capital_tiebreak)
                best=idx
                objective=score
                capital_tiebreak=distance
            end
        end
    end
    if best==0
        minx=minimum(root.X for root in enumeration.roots)
        maxx=maximum(root.X for root in enumeration.roots)
        throw(DomainError((
            t=t,J=J,n_fixed=length(enumeration.roots),
            min_X=minx,max_X=maxx,expected_X=expected.X_next,
            original_productivity_draw_not_rejected=true
        ),"No positive-investment multi-step equilibrium; dynamic path terminates"))
    end
    root=enumeration.roots[best]

    # Regression: original J=1 path through t=4 and the UNCHANGED t=5
    # productivity profile must reproduce Stage12's exact J=2 equilibrium.
    if J==2 && t==5 && state.base_seed==20261101 &&
       state.replication==1 && length(state.firms.capital)==110000
        isapprox(root.K,STAGE13_REFERENCE_K5;atol=2e-9,rtol=1e-10) ||
            error("original t5 Stage12 J=2 capital regression failed")
        isapprox(root.X,STAGE13_REFERENCE_X5;atol=2e-9,rtol=1e-10) ||
            error("original t5 Stage12 J=2 investment regression failed")
    end

    # The selected root's direct firm-level response is independently
    # checked, but the full set of roots is not reaggregated every period.
    direct=multistep_direct_response(state.firms,cal,root.K,J)
    isapprox(direct.K,root.K;atol=2e-9,rtol=1e-10) ||
        error("selected multi-step root failed independent capital check")
    isapprox(direct.X,root.X;atol=2e-9,rtol=1e-10) ||
        error("selected multi-step root failed independent investment check")

    firms=state.firms
    workspace=state.workspace
    n=length(firms.capital)
    total_X=0.0
    invest_one=0
    invest_many=0
    divest_one=0
    divest_many=0
    inactive=0
    @inbounds for i in 1:n
        d=(1-p.delta)*firms.capital[i]
        li=firms.lambda[i]
        b=threshold_coefficient(li,a_next[i],cal)
        action=multistep_policy_at_K(root.K,d,li,b,J)
        next_k=d*li^action
        workspace.next_capital[i]=next_k
        total_X+=next_k-d
        invest_one+=(action==1)
        invest_many+=(action>=2)
        divest_one+=(action==-1)
        divest_many+=(action<=-2)
        inactive+=(action==0)
    end
    invest_one+invest_many+divest_one+divest_many+inactive==n ||
        error("firm investment action counts do not sum to N")
    invest_many==direct.invest_many && divest_many==direct.divest_many ||
        error("firm actions disagree with independent policy evaluation")
    X_direct=total_X/n
    K_direct=aggregate_capital(a_next,workspace.next_capital,p.alpha,rho)
    Kres=K_direct-root.K
    Xres=X_direct-root.X
    abs(Kres)<=1e-9*max(1,abs(root.K)) ||
        error("exact next firm CES does not match selected capital root")
    abs(Xres)<=1e-9*max(1,abs(root.X)) ||
        error("exact next firm investment does not match root X")

    # All macro calculations follow the original first-order continuation,
    # not a newly solved nonlinear rational-expectations equilibrium.
    logK_prior=log(kp/steady.K)
    logK_next=log(root.K/steady.K)
    logmu_prior=log(mup/steady.mu)
    logmu=state.re.h*logK_next
    mu=steady.mu*exp(logmu)
    logX=log(root.X/steady.X)
    isfinite(mu) && mu>0 || error("invalid stage13 marginal-utility continuation")
    realized=Nirei2015._realized_macroeconomy(
        cal,state.re,logK_prior,logX,logmu,logmu_prior
    )
    eps=log(root.X/expected.X_next)
    logK_linear=expected.logK_next+p.delta*eps
    level_gap=(realized.Y-realized.C-root.X)/realized.Y
    row=Dict{String,Any}(
        "t"=>t,"J"=>J,"K_prior"=>kp,"K_next"=>root.K,
        "X"=>root.X,"Y"=>realized.Y,"C"=>realized.C,"L"=>realized.L,
        "mu"=>mu,"epsilon"=>eps,"expected_K"=>expected.K_next,
        "expected_X"=>expected.X_next,
        "n_fixed"=>length(enumeration.roots),
        "n_positive_X"=>positives,
        "invest_1"=>invest_one,"invest_2plus"=>invest_many,
        "divest_1"=>divest_one,"divest_2plus"=>divest_many,
        "inactive"=>inactive,"CES_residual"=>Kres,
        "investment_residual"=>Xres,
        "linear_K_gap"=>logK_next-logK_linear,
        "logK_forecast_error"=>logK_next-expected.logK_next,
        "logmu_forecast_error"=>logmu-expected.logmu_next,
        "stable_jump_residual"=>
            (logmu-expected.logmu_next)-state.re.h*(logK_next-expected.logK_next),
        "level_resource_gap"=>level_gap,
        "linear_resource_residual"=>realized.resources_residual,
        "linear_euler_residual"=>realized.euler_residual,
        "firm_count"=>n,"shock_checksum"=>shock_checksum
    )

    # Commit only when root and macro checks have passed.
    previous_capital=firms.capital
    firms.capital=workspace.next_capital
    workspace.next_capital=previous_capital
    advance_productivity!(firms)
    state.K=root.K
    state.mu_prior=mu
    state.period=t+1
    return row
end

function stage13_single_experiment(experiment::String,T::Int,N::Int,seed::Int)
    experiment in ("A","B") || error("Expected experiment A or B")
    state=initialize_simulation(;N=N,base_seed=seed)
    rows=Dict{String,Any}[]
    failed_t=-1
    fail_reason=""
    for t in 0:(T-1)
        try
            row=if experiment=="A" && t<5
                period=simulate_period!(state)
                stage13_period_row(period,N)
            else
                stage13_multistep_step!(state,2)
            end
            row["experiment"]=experiment
            push!(rows,row)
            println("STAGE13 ",experiment," t=",t,
                    " J=",row["J"]," K=",row["K_next"]," X=",row["X"],
                    " feasible_roots=",row["n_positive_X"])
        catch err
            err isa InterruptException && rethrow()
            failed_t=t
            fail_reason=sprint(showerror,err)
            println("STAGE13 ",experiment," STOPPED t=",t," reason=",fail_reason)
            break
        end
    end
    length(rows)==state.period ||
        error("stored rows do not match committed transition dates")
    status=length(rows)==T ? "completed" : "failed"
    summary=Dict{String,Any}(
        "experiment"=>experiment,"status"=>status,"target_periods"=>T,
        "completed_periods"=>length(rows),"first_failure_t"=>failed_t,
        "failure_reason"=>fail_reason,"N"=>N,"seed"=>seed,
        "initial_K"=>isempty(rows) ? state.K : first(rows)["K_prior"],
        "final_K"=>state.K,"final_mu"=>state.mu_prior,
        "max_CES_residual"=>isempty(rows) ? 0.0 :
            maximum(abs(row["CES_residual"]) for row in rows),
        "max_X_residual"=>isempty(rows) ? 0.0 :
            maximum(abs(row["investment_residual"]) for row in rows),
        "max_abs_level_resource_gap"=>isempty(rows) ? 0.0 :
            maximum(abs(row["level_resource_gap"]) for row in rows),
        "min_X"=>isempty(rows) ? "undefined" :
            minimum(row["X"] for row in rows),
        "min_K"=>isempty(rows) ? state.K :
            minimum(row["K_next"] for row in rows),
        "max_K"=>isempty(rows) ? state.K :
            maximum(row["K_next"] for row in rows),
        "mean_invest_2plus_share"=>isempty(rows) ? "undefined" :
            mean(row["invest_2plus"]/N for row in rows),
        "model"=>"J=2 diagnostic counterfactual with original linear continuum forecast",
        "is_original_paper_replication"=>false,
        "no_shock_rejection"=>true
    )
    if length(rows)==1100
        summary["mean_log_X"]=mean(log(row["X"]) for row in rows)
        summary["std_log_X"]=std(log(row["X"]) for row in rows)
        summary["mean_log_Y"]=mean(log(row["Y"]) for row in rows)
        summary["std_log_Y"]=std(log(row["Y"]) for row in rows)
    end
    return (rows=rows,summary=summary)
end

function stage13_output(experiment,T,N,seed,output)
    path=normpath(joinpath(@__DIR__,"..","results"))
    mkpath(path)
    stem="diagnostic_multistep_simulation_$(experiment)_T$(T)_N$(N)_seed$(seed)"
    cols=(
        "experiment","t","J","K_prior","K_next","X","Y","C","L","mu",
        "epsilon","expected_K","expected_X","n_fixed","n_positive_X",
        "invest_1","invest_2plus","divest_1","divest_2plus","inactive",
        "CES_residual","investment_residual","linear_K_gap",
        "logK_forecast_error","logmu_forecast_error","stable_jump_residual",
        "level_resource_gap","linear_resource_residual",
        "linear_euler_residual","firm_count","shock_checksum"
    )
    stage12_csv(joinpath(path,stem*".csv"),output.rows,cols)
    open(joinpath(path,stem*".toml"),"w") do io
        TOML.print(io,output.summary)
    end
    println("STAGE13 SUMMARY ",output.summary)
    println("STAGE13 FILES ",stem,".{csv,toml}")
end

function stage13_run(T::Int,N::Int,seed::Int)
    T>=6 && N>0 && seed>=0 ||
        error("T must be >= 6, N must be positive and seed nonnegative")
    allresults=Dict{String,Any}()
    for experiment in ("A","B")
        result=stage13_single_experiment(experiment,T,N,seed)
        stage13_output(experiment,T,N,seed,result)
        allresults[experiment]=result
    end
    # Period-5 raw productivity must be the same in A and B despite different
    # policies in their prior states. Check at t>=5 when both reached date 5.
    ar=allresults["A"].rows
    br=allresults["B"].rows
    overlap=min(length(ar),length(br))
    for j in 1:overlap
        t=j-1
        if t>=5
            ar[j]["shock_checksum"]==br[j]["shock_checksum"] ||
                error("A/B productivity RNG diverged at date t=$t")
        end
    end
    return allresults
end

function stage13_smoke()
    stage12_small_tests()
    # An initial J=2 step from a small seeded state checks capital and
    # productivity-buffer updating without invoking the costly N=110000 sweep.
    state=initialize_simulation(;N=48,base_seed=20261101)
    firstrow=stage13_multistep_step!(state,2)
    state.period==1 && !state.firms.next_ready ||
        error("small-N firm date/productivity buffers did not advance")
    isapprox(aggregate_capital(state.firms,state.calibration),state.K;
             atol=2e-9,rtol=1e-10) ||
        error("small-N aggregate K does not match committed firm arrays")
    firstrow["X"]>0.0 || error("small-N smoke has non-positive selected investment")
    println("STAGE13 SMALL DYNAMIC SMOKE PASS firstK=",state.K,
            " firstX=",firstrow["X"])
end

if abspath(PROGRAM_FILE)==@__FILE__
    mode=isempty(ARGS) ? "smoke" : ARGS[1]
    if mode=="smoke"
        stage13_smoke()
    elseif mode=="run"
        T=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 20
        N=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 110000
        seed=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 20261101
        stage13_run(T,N,seed)
    else
        error("Usage: diagnostic_multistep_simulation.jl smoke | run [T] [N] [seed]")
    end
end
