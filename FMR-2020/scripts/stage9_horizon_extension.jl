# Stage 9 — two-branch T=200 -> T=300 horizon convergence.
#
# Uses the existing Euler-polished strict-portfolio 3T-1 GE residual and solver,
# and an IDENTICAL original 200-date pre-shock equilibrium, shock and initial
# household distribution. No change to consumption, distribution or market
# equations. Source T=200 roots are read-only.
#
# Run from FMR-2020:
#  julia --project=. scripts/stage9_horizon_extension.jl --audit
#  julia --project=. scripts/stage9_horizon_extension.jl --solve-depreciation
#  julia --project=. scripts/stage9_horizon_extension.jl --solve-appreciation
#  julia --project=. scripts/stage9_horizon_extension.jl --compare
#
# Only --audit is permitted in CI. The 899-column nonlinear solver is used
# by explicit opt-in on the owner PC, without artificial 5-step limits.

using FMR2020
using LinearAlgebra

include(joinpath(@__DIR__,"stage8_compare_saved_roots.jl"))

const STAGE9_HORIZON=300
const STAGE9_SHOCK_DATE=41
const STAGE9_LEVERAGE=1.0/16.0

"""
    stage9_assert_same_economy(old, extended)

Ensure the ONLY model change is the terminal truncation date. Any change
in inherited gross portfolios, household distribution, grid, exogenous
credit path through date 200, productivity, or terminal steady state
is an error (not a horizon-sensitivity experiment).
"""
function stage9_assert_same_economy(old::Stage3Context{T},
                                     extended::Stage3Context{T}) where {T<:AbstractFloat}
    a=old.economy
    b=extended.economy
    n=length(a.B)
    m=length(b.B)
    m>n || throw(DimensionMismatch("extended horizon must grow"))
    a.B==b.B[1:n] ||
        throw(ArgumentError("foreign-credit shock changed at old dates"))
    a.z==b.z[1:n] ||
        throw(ArgumentError("aggregate productivity changed at old dates"))
    all(v->v==a.terminal.B,b.B[n+1:end]) ||
        throw(ArgumentError("additional foreign-credit supply is nonstationary"))
    all(v->v==a.z[end],b.z[n+1:end]) ||
        throw(ArgumentError("additional productivity is not stationary"))
    isequal(a.initial,b.initial) ||
        throw(ArgumentError("inherited initial economy changed"))
    isequal(a.terminal,b.terminal) ||
        throw(ArgumentError("terminal stationary equilibrium changed"))
    isequal(a.params,b.params) && isequal(a.chain,b.chain) &&
    isequal(a.grids,b.grids) && isequal(a.reference,b.reference) ||
        throw(ArgumentError("parameters, Markov chain, asset grids or reference changed"))
    old.mu_pre==extended.mu_pre ||
        throw(ArgumentError("inherited household joint distribution changed"))
    old.terms_pre==extended.terms_pre &&
    old.leverage_hat==extended.leverage_hat &&
    old.foreign==extended.foreign && old.theta_nr==extended.theta_nr ||
        throw(ArgumentError("portfolio/revaluation/pricing regime changed"))
    sum(old.mu_pre)>zero(T) &&
        abs(sum(old.mu_pre)-one(T))<=T(1e-8) ||
        throw(ArgumentError("inherited probability mass is invalid"))
    all(>=(zero(T)),extended.mu_pre) ||
        throw(ArgumentError("inherited distribution contains negative probabilities"))
    return (
        old_horizon=n,new_horizon=m,
        old_unknown_count=3n-1,new_unknown_count=3m-1,
        inherited_K=a.initial.K,
        inherited_B=a.B_initial,
        inherited_q=a.q_initial,
        terms_pre=old.terms_pre,
        inherited_distribution_mass=sum(extended.mu_pre),
        credit_path_prefix_max_difference=maximum(abs,b.B[1:n]-a.B),
        same_terminal=true,same_grid=true,same_initial_distribution=true,
    )
end

"""
    stage9_extend_saved_root(old_path, old_ctx, new_ctx)

Numerical initial guess ONLY. The original 200-period K/L/terms path
is copied without modification; only periods 201–300 are appended using
the existing canonical steady-state prolongation.
"""
function stage9_extend_saved_root(old_path::AbstractVector{T},
                                  old_ctx::Stage3Context{T},
                                  new_ctx::Stage3Context{T}) where {T<:AbstractFloat}
    stage9_assert_same_economy(old_ctx,new_ctx)
    n=length(old_ctx.economy.B)
    m=length(new_ctx.economy.B)
    length(old_path)==3n-1 ||
        throw(DimensionMismatch("old market-root length differs"))
    ans=FMR2020.stage3_prolong_aggregate_guess(
        old_path,old_ctx,new_ctx)
    length(ans)==3m-1 ||
        throw(DimensionMismatch("new market candidate must contain 899 variables"))
    ans[1:n-1]==old_path[1:n-1] ||
        throw(ErrorException("old capital path changed during prolongation"))
    ans[m: m+n-1]==old_path[n:2n-1] ||
        throw(ErrorException("old labor path changed during prolongation"))
    ans[2m:2m+n-1]==old_path[2n:3n-1] ||
        throw(ErrorException("old exchange-price path changed"))
    return ans
end

function stage9_filepaths(base::AbstractString)
    dir=joinpath(base,"results","stage9_horizon")
    prior=joinpath(base,"results","economic_reproduction")
    stage8=joinpath(base,"results","stage8_depreciation")
    oldroot=(
        appreciation=joinpath(prior,
            "economic_flexible_leveraged_euler_polished_"*
            "STRICT_PORTFOLIO_MARKET_ROOT_UNCERTIFIED.csv"),
        depreciation=joinpath(stage8,
            "source_initialized_STRICT_MARKET_ROOT_UNCERTIFIED.csv"),
    )
    return (
        dir=dir,prior=prior,oldroot=oldroot,
        newroot=(
            appreciation=joinpath(dir,
                "appreciation_T300_MARKET_ROOT_UNCERTIFIED.csv"),
            depreciation=joinpath(dir,
                "depreciation_T300_MARKET_ROOT_UNCERTIFIED.csv"),
        ),
    )
end

"""
    stage9_contexts(steady, preboom)

Reconstruct ONE strict inherited context from the *existing saved T=200
Euler-polished expansion market root*. It is not recalibrated for T=300.
"""
function stage9_contexts(steady,preboom)
    old=FMR2020.stage3_reversal_context(
        steady,preboom;shock_start=STAGE9_SHOCK_DATE,
        leverage_hat=STAGE9_LEVERAGE)
    length(old.economy.B)==200 ||
        throw(DimensionMismatch("Stage 8 comparison must start at T=200"))
    long=FMR2020.stage3_extend_horizon(old,STAGE9_HORIZON)
    inv=stage9_assert_same_economy(old,long)
    return old,long,inv
end

"""
    stage9_terminal_summary(work, context)

Independent terminal copied-policy budget diagnostic on every strictly
positive-mass state, without truncation, renormalization or state filtering.
The additional max on cells >1e-12 is informational only.
"""
function stage9_terminal_summary(w::FMR2020.TransitionWorkspace{T},
                                 ctx::Stage3Context{T}) where {T<:AbstractFloat}
    n=length(ctx.economy.B)
    ss=ctx.economy.terminal
    actual_m=FMR2020.consumption_cost(w.terms[n],ctx.economy.params)
    ss_m=FMR2020.consumption_cost(ss.terms,ctx.economy.params)
    grid=ctx.economy.grids.fine
    masses=T[]
    budgets=T[]
    worst_decomp=zero(T)
    worst_ss=zero(T)
    for j in axes(w.distribution,2),i in eachindex(grid)
        mass=w.distribution[i,j,n]
        mass>zero(T) || continue
        z=stage8_terminal_budget_components(
            grid[i],ctx.economy.chain.levels[j],
            w.fine_consumption[i,j,n],w.fine_savings[i,j,n],
            w.rate[n],w.wage[n]*w.L[n],w.profits[n],actual_m,
            ss.r,ss.wage*ss.L,ss.profits,ss_m)
        worst_decomp=max(worst_decomp,abs(z.decomposition_error))
        worst_ss=max(worst_ss,abs(z.stationary_policy_budget))
        push!(masses,mass)
        push!(budgets,z.actual_budget)
    end
    tail=stage8_terminal_mass_summary(masses,budgets)
    worst_decomp<=T(1e-9) ||
        throw(ErrorException("independent terminal budget identity does not hold"))
    return merge(tail,(
        stationary_policy_budget_max=worst_ss,
        decomposition_max_error=worst_decomp,
        terminal_K_gap=w.K[n]-ss.K,
        terminal_r_gap=w.rate[n]-ss.r,
        terminal_wage_gap=w.wage[n]-ss.wage,
        terminal_profits_gap=w.profits[n]-ss.profits,
        terminal_L_gap=w.L[n]-ss.L,
        terminal_terms_gap=w.terms[n]-ss.terms,
    ))
end

function stage9_save_metrics(path,values)
    stage8_store_metrics(path,values)
end

function stage9_audit_guess(ctx::Stage3Context{T},
                            x::Vector{T},name::String) where {T<:AbstractFloat}
    eval=stage8_one_residual(ctx,x;partials=false)
    rows=Pair{String,Any}[
        "branch"=>name,"horizon"=>length(ctx.economy.B),
        "unknowns"=>length(x),
        "pre_shock_terms"=>ctx.terms_pre,
        "impact_terms"=>x[2length(ctx.economy.B)],
        "nominal_depreciation_logpoints"=>
            100log(ctx.terms_pre/x[2length(ctx.economy.B)]),
        "candidate_status"=>(eval.ok ? "DEFINED_NOT_EQUILIBRIUM" :
                                        "ECONOMICALLY_UNDEFINED"),
    ]
    append!(rows,eval.metrics)
    if eval.ok
        for (k,v) in pairs(stage9_terminal_summary(eval.work,ctx))
            push!(rows,"terminal_"*string(k)=>v)
        end
    end
    return eval,rows
end

function stage9_bootstrap()
    base=normpath(joinpath(@__DIR__,".."))
    paths=stage9_filepaths(base)
    steady=FMR2020.solve_three_steady_states(
        FMR2020.default_parameters();verbose=false)
    preboom=stage5_existing_boom(steady,paths.prior)
    old,long,invariants=stage9_contexts(steady,preboom)
    oldpaths=(
        appreciation=read_economic_root_guess(paths.oldroot.appreciation,old),
        depreciation=read_economic_root_guess(paths.oldroot.depreciation,old),
    )
    oldpaths.appreciation===nothing &&
        error("Missing original Stage 8 appreciation market-root CSV")
    oldpaths.depreciation===nothing &&
        error("Missing original Stage 8 depreciation market-root CSV")
    n=length(old.economy.B)
    oldpaths.appreciation[2n]>old.terms_pre ||
        error("Saved appreciation root has the wrong FX sign")
    oldpaths.depreciation[2n]<old.terms_pre ||
        error("Saved depreciation root has the wrong FX sign")
    guesses=(
        appreciation=stage9_extend_saved_root(oldpaths.appreciation,old,long),
        depreciation=stage9_extend_saved_root(oldpaths.depreciation,old,long),
    )
    return (old=old,long=long,paths=paths,
            invariants=invariants,oldpaths=oldpaths,guesses=guesses)
end

function stage9_audit(b;save::Bool=true)
    rows=Pair{String,Any}[]
    for (k,v) in pairs(b.invariants)
        push!(rows,string(k)=>v)
    end
    push!(rows,"same_foreign_credit_shock"=>true)
    push!(rows,"same_pre_shock_joint_distribution"=>true)
    push!(rows,"same_terminal_stationary_equilibrium"=>true)
    push!(rows,"199_to_299_asset_equations"=>true)
    if save
        stage9_save_metrics(
            joinpath(b.paths.dir,"horizon_initialization_audit.csv"),rows)
    end
    summaries=Pair{String,Any}[]
    evaluations=Dict{String,Any}()
    for branch in ("depreciation","appreciation")
        guess=getproperty(b.guesses,Symbol(branch))
        e,report=stage9_audit_guess(b.long,guess,branch)
        evaluations[branch]=e
        if save
            stage9_save_metrics(
                joinpath(b.paths.dir,branch*"_T300_initial_residuals.csv"),
                report)
        end
        for (key,value) in report
            push!(summaries,branch*"_"*key=>value)
        end
        println("STAGE9 INITIAL ",branch,
                " terms=",guess[2STAGE9_HORIZON],
                " full899_defined=",e.ok,
                " max_scaled=",e.ok ? maximum(abs,e.residual) : "UNDEFINED")
    end
    if save
        stage9_save_metrics(
            joinpath(b.paths.dir,"horizon_initial_residuals.csv"),summaries)
    end
    return evaluations
end

function stage9_solve_one(b,branch::String)
    branch in ("appreciation","depreciation") ||
        throw(ArgumentError("Stage 9 requested unknown branch"))
    n=length(b.long.economy.B)
    n==STAGE9_HORIZON || throw(DimensionMismatch("expected T=300"))
    rootfile=getproperty(b.paths.newroot,Symbol(branch))
    checkpoint=joinpath(b.paths.dir,branch*"_T300_UNCONVERGED.checkpoint")
    trace=joinpath(b.paths.dir,branch*"_T300_accepted_trace.csv")
    if isfile(rootfile)
        x=read_economic_root_guess(rootfile,b.long)
        x===nothing && error("Unreadable existing root")
        root_x1=x[2n]
        correct=branch=="depreciation" ?
            root_x1<b.long.terms_pre :
            root_x1>b.long.terms_pre
        correct || error("Saved Stage 9 root belongs to the opposite FX branch")
        e,rows=stage9_audit_guess(b.long,x,branch)
        e.ok && maximum(abs,e.residual)<=1e-6 ||
            error("Previously saved Stage 9 market root does not solve current equations")
        stage9_save_metrics(
            joinpath(b.paths.dir,branch*"_T300_verified_existing.csv"),rows)
        println("STAGE9 VERIFIED EXISTING ",branch,
                " market root; no nonlinear GE solve or overwrite")
        return nothing
    end
    start=getproperty(b.guesses,Symbol(branch))
    previous=0
    if isfile(checkpoint)
        # Incomplete iterations are never confused with converged roots.
        saved=stage8_read_checkpoint(checkpoint,3n-1)
        e=stage8_one_residual(b.long,saved.x;partials=false)
        e.ok || error("Existing Stage 9 checkpoint infeasible: "*string(e.metrics))
        stage8_verify_resume_residual(
            saved.stored_max,maximum(abs,e.residual))
        start=saved.x
        previous=saved.iterations
        println("STAGE9 VERIFIED RESUME ",branch,
                " total accepted=",previous,
                " max scaled=",saved.stored_max)
    else
        e=stage8_one_residual(b.long,start;partials=false)
        e.ok || error("Stage 9 initial guess infeasible: "*string(e.metrics))
    end
    # Snapshot at each solver-accepted step; never change the Stage 8 files.
    function accepted_step(x,f,i)
        total=previous+i
        snapshot=joinpath(b.paths.dir,branch*"_T300_iteration_"*
                          lpad(string(total),4,'0')*".checkpoint")
        !isfile(snapshot) ||
            error("Refusing to overwrite accepted iterate "*snapshot)
        FMR2020.stage3_write_checkpoint(snapshot,x,f,total)
        FMR2020.stage3_write_checkpoint(checkpoint,x,f,total)
        stage8_append_accepted_step(
            trace,"Stage9-"*branch,previous,i,
            x,f,b.long.terms_pre,n)
        println("STAGE9 ACCEPTED ",branch,
                " total=",total," impact_terms=",x[2n],
                " depreciation=",x[2n]<b.long.terms_pre,
                " max899=",maximum(abs,f))
        return nothing
    end
    println("STAGE9 EXISTING SOLVER T300 ",branch,
            " maxiter=70 max_evaluations=5000; no artificial 5-step cutoff")
    sol=FMR2020.solve_stage3_flexible_full(
        b.long;x0=copy(start),strict_portfolio=true,
        polish_households=true,tolerance=1e-6,
        maxiter=70,max_evaluations=5000,
        on_accepted_step=accepted_step,verbose=true)
    sol.converged || error("Stage 9 did not converge")
    # Independent current-equation and household/economic checks.
    work=FMR2020.TransitionWorkspace(b.long.economy)
    residual=zeros(eltype(sol.x),3n-1)
    FMR2020.stage3_flexible_full_residual!(
        residual,sol.x,b.long,work;
        strict_portfolio=true,polish_households=true)
    maximum(abs,residual)<=1e-6 ||
        error("Independently re-evaluated 899 root fails market clearing")
    correct_branch=branch=="depreciation" ?
        work.terms[1]<b.long.terms_pre :
        work.terms[1]>b.long.terms_pre
    correct_branch || error(
        "899 solver reached a valid candidate on the OPPOSITE FX branch; "*
        "preserving accepted checkpoints, but refusing to mislabel "*
        branch*" as a separate "*branch*" market root. x1="*
        string(work.terms[1]))
    checked=FMR2020.TransitionSolveResult(
        b.long.economy,work,sol.x,residual,0,
        work.function_evals,maximum(abs,residual),true)
    diag=FMR2020.stage3_diagnostics(
        checked,b.long;explicit_flexible=true)
    term=stage9_terminal_summary(work,b.long)
    market=(diag.max_asset_market_level<=1e-4 &&
            diag.max_labor_market_level<=1e-4 &&
            diag.max_goods_market_level<=1e-4 &&
            diag.max_household_budget<=1e-8 &&
            diag.max_mass_deviation<=1e-8)
    market || error("Returned Stage 9 root fails independent market / budget tests")
    terminal_certified=(term.max_budget_all_occupied<=1e-6 &&
                        abs(term.terminal_K_gap)/max(1,abs(b.long.economy.terminal.K))<=1e-4)
    metrics=Pair{String,Any}[
        "branch"=>branch,"horizon"=>n,
        "max_899_scaled_residual"=>maximum(abs,residual),
        "impact_terms"=>work.terms[1],
        "nominal_depreciation_logpoints"=>
            100log(b.long.terms_pre/work.terms[1]),
        "finite_horizon_market_root_verified"=>true,
        "infinite_horizon_certified"=>false,
        "terminal_household_certified"=>terminal_certified,
        "max_euler_KKT"=>diag.max_dynamic_euler,
        "max_budget_interior"=>diag.max_household_budget,
        "max_mass_error"=>diag.max_mass_deviation,
        "max_asset_market_level"=>diag.max_asset_market_level,
        "max_labor_market_level"=>diag.max_labor_market_level,
        "max_goods_market_level"=>diag.max_goods_market_level,
    ]
    for (k,v) in pairs(term)
        push!(metrics,"terminal_"*string(k)=>v)
    end
    stage9_save_metrics(
        joinpath(b.paths.dir,branch*"_T300_economic_certification.csv"),
        metrics)
    # Use an explicit uncertified filename and never overwrite old 200 roots.
    FMR2020.stage3_write_csv(checked,b.long,rootfile)
    println("STAGE9 VERIFIED NEW FINITE-HORIZON MARKET ROOT ",branch,
            " saved=",rootfile," max899=",maximum(abs,residual),
            " terminal_budget=",term.max_budget_all_occupied,
            " impact_terms=",work.terms[1])
    return nothing
end

function stage9_compare(b)
    n=length(b.long.economy.B)
    n0=length(b.old.economy.B)
    results=Pair{String,Any}[
        "same_inherited_economy"=>true,
        "T200"=>n0,"T300"=>n,
        "infinite_horizon_multiplicity_certified"=>false,
    ]
    impacts=Dict{String,Float64}()
    for branch in ("depreciation","appreciation")
        before=getproperty(b.oldpaths,Symbol(branch))
        file=getproperty(b.paths.newroot,Symbol(branch))
        after=read_economic_root_guess(file,b.long)
        after===nothing && error("Missing Stage 9 300-period root for "*branch)
        e=stage8_one_residual(b.long,after;partials=false)
        e.ok && maximum(abs,e.residual)<=1e-6 ||
            error("T300 "*branch*" is not a full 899 market root")
        oldn=n0
        # Recalculate old and new policy & distribution using one common
        # model definition at each horizon, not just stored CSV consumption.
        eold=stage8_one_residual(b.old,before;partials=false)
        eold.ok && maximum(abs,eold.residual)<=1e-6 ||
            error("T200 "*branch*" is no longer a current GE root")
        w1=eold.work
        w2=e.work
        term200=stage9_terminal_summary(w1,b.old)
        term300=stage9_terminal_summary(w2,b.long)
        impact200=w1.terms[1]
        impact300=w2.terms[1]
        impacts[branch]=impact300
        push!(results,branch*"_x1_T200"=>impact200)
        push!(results,branch*"_x1_T300"=>impact300)
        push!(results,branch*"_impact_terms_change"=>impact300-impact200)
        push!(results,branch*"_nominal_depreciation_T200"=>
              100log(b.old.terms_pre/impact200))
        push!(results,branch*"_nominal_depreciation_T300"=>
              100log(b.long.terms_pre/impact300))
        push!(results,branch*"_full899_max_scaled"=>maximum(abs,e.residual))
        for horizon in (200,300)
            t=(horizon==200 ? term200 : term300)
            for (key,val) in pairs(t)
                push!(results,branch*"_T"*string(horizon)*"_"*string(key)=>val)
            end
        end
        # Only the common CALENDAR dates are compared.
        for horizon in (40,100,199)
            for (key,values1,values2) in (
                ("K",w1.K,w2.K),("L",w1.L,w2.L),
                ("terms",w1.terms,w2.terms),
                ("Y",w1.Y,w2.Y),
                ("investment",w1.I,w2.I),
                ("home_consumption",w1.aggregate_home,w2.aggregate_home),
            )
                errorval=maximum(
                    abs(values2[i]-values1[i])/
                    max(1,abs(values1[i])) for i in 1:horizon)
                push!(results,branch*"_"*key*"_common_first"*
                      string(horizon)*"_max_relative_gap"=>errorval)
            end
        end
        println("STAGE9 HORIZON ROOT ",branch,
                " x1 T200=",impact200,
                " x1 T300=",impact300,
                " terminal weighted abs T200=",
                term200.weighted_absolute_budget,
                " T300=",term300.weighted_absolute_budget,
                " material max T200=",term200.max_budget_material_cells,
                " T300=",term300.max_budget_material_cells)
    end
    impacts["depreciation"]<b.long.terms_pre<
        impacts["appreciation"] ||
        error("Two distinct opposite exchange-rate branches did not survive")
    push!(results,"two_distinct_T300_finite_market_roots"=>true)
    stage9_save_metrics(
        joinpath(b.paths.dir,"horizon_convergence_comparison.csv"),
        results)
    println("STAGE9 BOTH 899 MARKET ROOTS IN SAME MODEL: verified; ",
            "terminal and infinite-horizon certification still separate")
    return results
end

function stage9_main()
    mode=isempty(ARGS) ? "--audit" : only(ARGS)
    mode in ("--audit","--solve-depreciation","--solve-appreciation",
             "--compare") ||
        error("Usage: julia --project=. scripts/stage9_horizon_extension.jl [--audit|--solve-depreciation|--solve-appreciation|--compare]")
    b=stage9_bootstrap()
    if mode=="--audit"
        stage9_audit(b)
        println("STAGE9 AUDIT COMPLETE: both 899 starting candidates tested. No GE root solve.")
    elseif mode=="--compare"
        stage9_compare(b)
    else
        branch=replace(mode,"--solve-"=>"")
        stage9_solve_one(b,branch)
    end
    return nothing
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    stage9_main()
end
