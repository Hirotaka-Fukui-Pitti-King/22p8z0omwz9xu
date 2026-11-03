# Stage 8 independent verification of TWO saved strict-portfolio 599-market roots.
#
# Evaluates the same archived Julia inherited household distribution and
# foreign-credit path with BOTH saved aggregate paths. Does NOT search for
# a new GE root, alter either CSV, or claim infinite-horizon certification.
#
# Run from FMR-2020:
#   julia --project=. scripts/stage8_compare_saved_roots.jl
#
# Requires BOTH local untracked market-root CSVs from actual PC solves.

using FMR2020
using LinearAlgebra

include(joinpath(@__DIR__, "stage8_depreciation_continuation.jl"))

"""
    stage8_terminal_budget_components(asset, productivity, consumption, savings,
                                      rate_actual, labor_income_actual, profits_actual,
                                      cost_actual, rate_ss, labor_income_ss,
                                      profits_ss, cost_ss)

Exact algebraic difference between an actual terminal-period individual
budget and the stationary continuation budget, at the SAME terminal
wealth and stationary policy. It decomposes residuals without asserting
that the SS household budget necessarily holds at an off-grid wealth node.
"""
function stage8_terminal_budget_components(a::T, productivity::T,
        c::T, savings::T, r_actual::T, wage_labor_actual::T,
        profits_actual::T, m_actual::T, r_ss::T,
        wage_labor_ss::T, profits_ss::T, m_ss::T) where {T<:AbstractFloat}
    all(isfinite, (a,productivity,c,savings,r_actual,wage_labor_actual,
                   profits_actual,m_actual,r_ss,wage_labor_ss,profits_ss,m_ss)) ||
        throw(ArgumentError("Nonfinite terminal comparison input"))
    stationary=(one(T)+r_ss)*a +
               wage_labor_ss*productivity+profits_ss-m_ss*c-savings
    interest=(r_actual-r_ss)*a
    wage=(wage_labor_actual-wage_labor_ss)*productivity
    profits=profits_actual-profits_ss
    home_cost=-(m_actual-m_ss)*c
    observed=(one(T)+r_actual)*a +
             wage_labor_actual*productivity+profits_actual-m_actual*c-savings
    decomposition=stationary+interest+wage+profits+home_cost
    return (actual_budget=observed,stationary_policy_budget=stationary,
            interest_contribution=interest,wage_contribution=wage,
            profit_contribution=profits,home_cost_contribution=home_cost,
            decomposition_error=observed-decomposition)
end

"""
    stage8_audit_saved_market_root(ctx, saved_path)

Read only the saved aggregate path, but independently re-solve household
EGM (with Euler polishing) and forward distribution and recompute ALL 599
strict-portfolio market equations at that path. Construct a diagnostic
result *only after* independently verifying the market residual <=1e-6.
This does not invoke a nonlinear 599-dimensional root solve.
"""
function stage8_audit_saved_market_root(ctx::Stage3Context{T},
                                        saved_path::AbstractString) where {T<:AbstractFloat}
    guess=read_economic_root_guess(saved_path,ctx)
    guess===nothing && throw(ArgumentError(
        "Missing saved owner-PC market-root CSV: "*saved_path))
    n=length(ctx.economy.B)
    length(guess)==3n-1 || throw(DimensionMismatch("Expected 599-root CSV"))
    w=FMR2020.TransitionWorkspace(ctx.economy)
    res=zeros(T,length(guess))
    FMR2020.stage3_flexible_full_residual!(
        res,guess,ctx,w;polish_households=true,strict_portfolio=true)
    maxres=maximum(abs,res)
    maxres<=T(1e-6) || throw(ErrorException(
        "Saved candidate does NOT satisfy the CURRENT strict-portfolio "*
        "599 full-market system; maximum residual=$(maxres): "*saved_path))
    # The diagnostic function's converged flag is used ONLY after
    # the external residual check; this is NOT a fresh numerical solve.
    result=FMR2020.TransitionSolveResult(
        ctx.economy,w,guess,copy(res),0,w.function_evals,maxres,true)
    report=FMR2020.stage3_diagnostics(
        result,ctx;explicit_flexible=true)
    support=stage8_portfolio_support(w,ctx,w.terms[1])
    support.violating_mass==0 && support.offgrid_mass==0 ||
        throw(ErrorException("Positive-mass original states fail feasibility"))

    fine=ctx.economy.grids.fine
    mu=view(w.distribution,:,:,n)
    stationary=ctx.economy.terminal
    p=ctx.economy.params
    actual_m=FMR2020.consumption_cost(w.terms[end],p)
    ss_m=FMR2020.consumption_cost(stationary.terms,p)
    worst=(abs=-one(T),state=(0,0),mass=zero(T),
           components=nothing)
    max_decomp=zero(T)
    max_ss_policy_budget=zero(T)
    sum_ss_policy_budget=zero(T)
    sum_interest=zero(T)
    sum_wage=zero(T)
    sum_profit=zero(T)
    sum_home_cost=zero(T)
    sum_budget=zero(T)
    sum_abs_budget=zero(T)
    mass_total=zero(T)
    for j in axes(mu,2), i in axes(mu,1)
        mass=mu[i,j]
        mass>zero(T) || continue
        c=w.fine_consumption[i,j,n]
        save=w.fine_savings[i,j,n]
        parts=stage8_terminal_budget_components(
            fine[i],ctx.economy.chain.levels[j],c,save,
            w.rate[n],w.wage[n]*w.L[n],w.profits[n],actual_m,
            stationary.r,stationary.wage*stationary.L,
            stationary.profits,ss_m)
        max_decomp=max(max_decomp,abs(parts.decomposition_error))
        max_ss_policy_budget=max(max_ss_policy_budget,
                                 abs(parts.stationary_policy_budget))
        sum_ss_policy_budget+=mass*parts.stationary_policy_budget
        sum_interest+=mass*parts.interest_contribution
        sum_wage+=mass*parts.wage_contribution
        sum_profit+=mass*parts.profit_contribution
        sum_home_cost+=mass*parts.home_cost_contribution
        sum_budget+=mass*parts.actual_budget
        sum_abs_budget+=mass*abs(parts.actual_budget)
        mass_total+=mass
        if abs(parts.actual_budget)>worst.abs
            worst=(abs=abs(parts.actual_budget),state=(i,j),
                   mass=mass,components=parts)
        end
    end
    abs(mass_total-one(T))<=T(1e-8) ||
        throw(ErrorException("Terminal probability mass not conserved"))
    max_decomp<=T(1e-9) ||
        throw(ErrorException("Terminal budget decomposition failed"))
    abs(sum_abs_budget/mass_total-report.weighted_terminal_household_budget)<=T(1e-9) ||
        throw(ErrorException("Independent terminal budget differs from diagnostics"))
    abs(worst.abs-report.max_terminal_household_budget)<=T(1e-9) ||
        throw(ErrorException("Independent maximum terminal budget differs from diagnostics"))
    return (path=saved_path,x=guess,max_scaled_residual=maxres,
            report=report,support=support,
            worst_terminal_budget_state=worst.state,
            worst_terminal_budget_mass=worst.mass,
            worst_terminal_budget_components=worst.components,
            terminal_stationary_policy_budget_max=max_ss_policy_budget,
            terminal_budget_decomposition_error=max_decomp,
            terminal_budget_signed_mass=sum_budget,
            terminal_budget_weighted_abs=sum_abs_budget,
            terminal_stationary_policy_signed_mass=sum_ss_policy_budget,
            terminal_interest_signed_mass=sum_interest,
            terminal_wage_signed_mass=sum_wage,
            terminal_profits_signed_mass=sum_profit,
            terminal_home_cost_signed_mass=sum_home_cost)
end

function stage8_two_root_report(ctx::Stage3Context, appreciation, depreciation,
                                target_dir::AbstractString)
    n=length(ctx.economy.B)
    a=appreciation.x[2n]
    d=depreciation.x[2n]
    xpre=ctx.terms_pre
    a>xpre && d<xpre ||
        throw(ErrorException(
            "Both independently verified market roots must lie on "*
            "opposite sides of the SAME inherited pre-shock terms of trade"))
    distance=norm(appreciation.x-depreciation.x)
    distance>0 || throw(ErrorException("Identical aggregate roots"))
    metrics=Pair{String,Any}[
        "source"=>"two owner-PC saved market roots; independently recomputed strict 599 map",
        "same_inherited_context"=>true,
        "same_foreign_credit_path"=>true,
        "same_Euler_polished_household_operator"=>true,
        "same_strict_gross_portfolio_convention"=>true,
        "same_terminal_continuation"=>true,
        "horizon"=>n,
        "pre_shock_terms"=>xpre,
        "appreciation_x1"=>a,
        "depreciation_x1"=>d,
        "appreciation_log_depreciation_points"=>100*log(xpre/a),
        "depreciation_log_depreciation_points"=>100*log(xpre/d),
        "full_599_path_euclidean_separation"=>distance,
        "two_distinct_finite_horizon_market_roots"=>true,
        "infinite_horizon_multiple_equilibria_certified"=>false,
    ]
    for (label,z) in (("appreciation",appreciation),
                      ("depreciation",depreciation))
        push!(metrics,label*"_strict_market_scaled_max"=>
              z.max_scaled_residual)
        for (name,value) in pairs(z.report)
            if value isa Number
                push!(metrics,label*"_"*string(name)=>value)
            end
        end
        for (name,value) in pairs(z.support)
            push!(metrics,label*"_support_"*string(name)=>value)
        end
        push!(metrics,label*"_worst_terminal_asset_index"=>
              z.worst_terminal_budget_state[1])
        push!(metrics,label*"_worst_terminal_productivity_index"=>
              z.worst_terminal_budget_state[2])
        push!(metrics,label*"_worst_terminal_state_mass"=>
              z.worst_terminal_budget_mass)
        for (name,value) in pairs(z.worst_terminal_budget_components)
            push!(metrics,label*"_worst_terminal_"*string(name)=>value)
        end
        for key in (
            "terminal_stationary_policy_budget_max",
            "terminal_budget_decomposition_error",
            "terminal_budget_signed_mass",
            "terminal_budget_weighted_abs",
            "terminal_stationary_policy_signed_mass",
            "terminal_interest_signed_mass",
            "terminal_wage_signed_mass",
            "terminal_profits_signed_mass",
            "terminal_home_cost_signed_mass",
        )
            push!(metrics,label*"_"*key=>getproperty(z,Symbol(key)))
        end
        # Both are only FINITE-HORIZON market roots until the terminal
        # continuation audit passes. A small 599 residual is insufficient.
        push!(metrics,label*"_terminal_certified"=>
              z.report.max_terminal_household_budget<=1e-6)
    end
    out=joinpath(target_dir,"two_market_root_comparison.csv")
    stage8_store_metrics(out,metrics)
    println("STAGE8 TWO-ROOT INDEPENDENT MARKET CHECK COMPLETE")
    println("Appreciation x1=",a,
            " max 599 residual=",appreciation.max_scaled_residual,
            " max terminal budget=",appreciation.report.max_terminal_household_budget)
    println("Depreciation x1=",d,
            " max 599 residual=",depreciation.max_scaled_residual,
            " max terminal budget=",depreciation.report.max_terminal_household_budget)
    println("Same inherited economy: YES, reconstructed at both saved roots")
    println("Infinite-horizon multiple-equilibrium certification: NO")
    println("Wrote ",out)
    return metrics
end

function stage8_compare_saved_roots_main()
    isempty(ARGS) || error("Usage: julia --project=. scripts/stage8_compare_saved_roots.jl")
    base=normpath(joinpath(@__DIR__,".."))
    roots=joinpath(base,"results","economic_reproduction")
    output=joinpath(base,"results","stage8_depreciation")
    lowpath=joinpath(output,
       "source_initialized_STRICT_MARKET_ROOT_UNCERTIFIED.csv")
    highpath=joinpath(roots,
       "economic_flexible_leveraged_euler_polished_"*
       "STRICT_PORTFOLIO_MARKET_ROOT_UNCERTIFIED.csv")
    isfile(lowpath) && isfile(highpath) ||
        error("Both saved PC market-root CSVs are required. "*
              "No root-solving fallback or fabricated source paths.")
    steady=FMR2020.solve_three_steady_states(
        FMR2020.default_parameters();verbose=false)
    boom=stage5_existing_boom(steady,roots)
    ctx=FMR2020.stage3_reversal_context(
        steady,boom;shock_start=41,leverage_hat=1.0/16.0)
    high=stage8_audit_saved_market_root(ctx,highpath)
    low=stage8_audit_saved_market_root(ctx,lowpath)
    stage8_two_root_report(ctx,high,low,output)
    return nothing
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    stage8_compare_saved_roots_main()
end
