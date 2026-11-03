# Stage 5: first-period foreign balance-sheet accounting and scalar
# trade-price branches at the SAVED 599-equation equilibrium paths.
# No full non-linear GE solve, no MATLAB source numerical dependency.
#
# julia --project=. scripts/stage5_initial_fx_audit.jl
using FMR2020

include(joinpath(@__DIR__,"run_economic_replication.jl"))

# Direct support test BEFORE interpolation of post-shock net wealth onto
# the fine asset grid. This distinguishes model-level infeasibility at
# occupied original states from discretization-created infeasible nodes.
# Factor prices are held fixed: candidate scalar root, NOT full GE.
function stage5_exact_original_state_support(w,ctx,terms;
                                             strict_portfolio::Bool=false)
    inp=ctx.economy
    grid=inp.grids.fine
    mu=ctx.mu_pre
    lev_positive=zero(eltype(grid))
    gross_wealth=zero(eltype(grid))
    for j in axes(mu,2),i in eachindex(grid)
        m=mu[i,j]
        m==0 && continue
        if grid[i]>0
            lev_positive+=m*grid[i]
        end
    end
    lev=inp.B_initial==0 ? one(eltype(grid)) :
        (1-ctx.leverage_hat)*inp.initial.K*inp.q_initial/lev_positive
    kfixed=ctx.leverage_hat*inp.initial.K
    for j in axes(mu,2),i in eachindex(grid)
        m=mu[i,j]
        m==0 && continue
        gross_wealth+=m*max(lev*grid[i]+kfixed,zero(eltype(grid)))
    end
    scale=strict_portfolio ?
          inp.q_initial*inp.initial.K/gross_wealth : one(eltype(grid))
    ratio=ctx.terms_pre/terms
    R=1+w.rate[1]
    minresources=Inf
    minstate=(0,0)
    minpost=Inf
    total_violating_mass=0.0
    count_violating=0
    for j in axes(mu,2),i in eachindex(grid)
        mass=mu[i,j]
        mass==0 && continue
        a=grid[i]
        capital=scale*max(lev*a+kfixed,zero(eltype(grid)))
        bond=a-capital
        a_post=capital+bond*ratio
        income=w.wage[1]*w.L[1]*inp.chain.levels[j]+w.profits[1]
        resources=R*a_post+income-inp.grids.borrowing_limit
        if resources<minresources
            minresources=resources
            minstate=(i,j)
            minpost=a_post
        end
        if resources<=0
            total_violating_mass+=mass
            count_violating+=1
        end
    end
    return (minimum_resources=minresources,worst_original_state=minstate,
            worst_postshock_wealth=minpost,
            violating_original_mass=total_violating_mass,
            violating_original_states=count_violating,
            feasible=minresources>0)
end

function stage5_existing_boom(steady,root_dir::AbstractString)
    ctx=stage3_expansion_context(steady)
    csv=joinpath(root_dir,
        "economic_credit_expansion_euler_polished_MARKET_ROOT_UNCERTIFIED.csv")
    x=read_economic_root_guess(csv,ctx)
    x===nothing && error("Missing corrected 200-date boom: "*csv)
    n=length(ctx.economy.B)
    w=FMR2020.TransitionWorkspace(ctx.economy)
    res=zeros(eltype(x),3n-1)
    FMR2020.stage3_flexible_full_residual!(
        res,x,ctx,w;polish_households=true)
    maxres=maximum(abs,res)
    maxres<=1e-6 || error("Saved corrected boom DOES NOT clear current markets: "*string(maxres))
    result=FMR2020.TransitionSolveResult(
        ctx.economy,w,x,res,0,w.function_evals,maxres,true)
    return result
end

function stage5_first_period_case(
    steady,boom,root_dir::AbstractString,label::String;
    leverage_hat::Float64=1/16,fixed::Bool=false)

    ctx=stage3_reversal_context(steady,boom;
                               shock_start=41,leverage_hat=leverage_hat)
    name="economic_"*label*"_euler_polished_MARKET_ROOT_UNCERTIFIED.csv"
    x=read_economic_root_guess(joinpath(root_dir,name),ctx)
    x===nothing && error("Missing saved corrected surprise GE: "*name)
    inp=ctx.economy
    p=inp.params
    w=FMR2020.TransitionWorkspace(inp)
    if fixed
        FMR2020.stage3_fixed_prices!(w,ctx,x)
    else
        FMR2020.stage3_flexible_full_prices!(w,ctx,x)
    end
    n=length(inp.B)
    res=zeros(Float64,3n-1)
    FMR2020.stage3_household_backward!(w,ctx)
    FMR2020.stage3_euler_polish_households!(w,ctx)
    FMR2020.stage3_distribution_forward!(w,ctx)
    if fixed
        FMR2020.stage3_market_residual!(res,w,ctx;fixed=true)
    else
        FMR2020.stage3_flexible_full_market_residual!(res,w,ctx)
    end
    maxres=maximum(abs,res)
    maxres<=1e-6 || error("Saved surprise GE no longer clears markets: "*name)

    distribution=similar(ctx.mu_pre)
    info=FMR2020.revalue_foreign_debt!(
        distribution,ctx.mu_pre,inp.grids.fine,
        inp.initial.K,inp.q_initial,inp.B_initial,
        ctx.terms_pre,w.terms[1],leverage_hat;foreign=true)
    strict_distribution=similar(ctx.mu_pre)
    corrected=FMR2020.revalue_foreign_debt!(
        strict_distribution,ctx.mu_pre,inp.grids.fine,
        inp.initial.K,inp.q_initial,inp.B_initial,
        ctx.terms_pre,w.terms[1],leverage_hat;
        foreign=true,strict_portfolio=true)

    qK=inp.q_initial*inp.initial.K
    pre_target=qK+inp.B_initial
    implied_wealth=qK+inp.B_initial*ctx.terms_pre/w.terms[1]
    gross=1+w.rate[1]
    predicted_consumption_gap=p.chi*gross*(info.net_wealth_after-implied_wealth)
    actual_consumption_gap=w.aggregate_home[1]-w.home_budget[1]
    # This is a check of the exact household-budget aggregation identity,
    # not a new equilibrium equation or an arbitrary numerical threshold.
    isapprox(predicted_consumption_gap,actual_consumption_gap;
             atol=5e-5,rtol=1e-2) ||
        error("Initial FX portfolio accounting does not explain aggregate consumption gap: "*name)

    t=1
    expgoods=w.Y[t]-w.I[t]-w.adjustment[t]
    inflation=w.terms[t]/ctx.terms_pre
    nominal_cost=fixed ? ctx.theta_nr/2*(inflation-1)^2 : 0.0
    outwealth=inp.B[t]+w.q[t]*w.K[t+1]
    other=w.wage[t]*w.L[t]+w.profits[t]
    C=expgoods-nominal_cost-
        p.chi*(gross*qK-outwealth+other)
    E=-p.chi*gross*inp.B_initial*ctx.terms_pre
    D=inp.reference.export_shift
    roots=FMR2020.foreign_debt_trade_roots(C,E,D,p.theta_star)
    root_fn(terms)=C+E/terms-D/terms^p.theta_star

    function check_branch(value)
        try
            a=FMR2020.initial_revaluation_budget_support!(w,ctx,value)
            return (root=value,feasible=a.feasible,
                    minimum_resources=a.minresource,
                    worst_grid_state=a.minstate,
                    mass=a.occupied_mass,error="")
        catch err
            if err isa DomainError
                return (root=value,feasible=false,
                        minimum_resources=NaN,
                        worst_grid_state=(0,0),
                        mass=NaN,error=sprint(showerror,err))
            end
            rethrow(err)
        end
    end
    lower=check_branch(roots.low)
    upper=isnothing(roots.high) ? nothing : check_branch(roots.high)
    observed=check_branch(w.terms[1])
    # Test the counterfactual corrected portfolio on the SAME prices and
    # household policies: this is one market RESIDUAL EVALUATION, not GE.
    # The strict mode revalues positions that add exactly to qK and B.
    FMR2020.stage3_distribution_forward!(
        w,ctx;strict_portfolio=true)
    strict_res=zeros(Float64,3n-1)
    if fixed
        FMR2020.stage3_market_residual!(
            strict_res,w,ctx;fixed=true)
    else
        FMR2020.stage3_flexible_full_market_residual!(
            strict_res,w,ctx)
    end
    strict_home_budget_gap=w.aggregate_home[1]-w.home_budget[1]
    strict_ge_residual=maximum(abs,strict_res)
    # Check each conditional trade-price branch under both portfolio
    # allocations, without re-solving 599 aggregate equations.
    function check_strict_branch(value)
        try
            v=FMR2020.initial_revaluation_budget_support!(
                w,ctx,value;strict_portfolio=true)
            return (feasible=v.feasible,minresource=v.minresource,
                    error="")
        catch err
            if err isa DomainError
                return (feasible=false,minresource=NaN,
                        error=sprint(showerror,err))
            end
            rethrow(err)
        end
    end
    strict_lower=check_strict_branch(roots.low)
    strict_upper=isnothing(roots.high) ? nothing :
        check_strict_branch(roots.high)

    direct_lower=stage5_exact_original_state_support(
        w,ctx,roots.low)
    direct_strict_lower=stage5_exact_original_state_support(
        w,ctx,roots.low;strict_portfolio=true)
    direct_upper=stage5_exact_original_state_support(
        w,ctx,roots.high)
    direct_strict_upper=stage5_exact_original_state_support(
        w,ctx,roots.high;strict_portfolio=true)

    return (
        economy=label,maximum_scaled_residual=maxres,
        terms_pre=ctx.terms_pre,terms_observed=w.terms[1],
        observed_nominal_depreciation=100*log(ctx.terms_pre/w.terms[1]),
        qK=qK,actual_gross_capital=info.capital,
        gross_capital_gap=info.capital-qK,
        B_initial=inp.B_initial,actual_foreign_bonds=info.bonds,
        foreign_bond_gap=info.bonds-inp.B_initial,
        pre_distribution_wealth=info.net_wealth_before,
        pre_aggregate_target_wealth=pre_target,
        revalued_wealth=info.net_wealth_after,
        predicted_wealth=implied_wealth,
        predicted_consumption_gap=predicted_consumption_gap,
        actual_consumption_gap=actual_consumption_gap,
        strictly_cleared_capital=corrected.capital,
        strictly_cleared_bonds=corrected.bonds,
        strict_portfolio_scale=corrected.capital_normalization,
        strict_portfolio_capital_gap=corrected.capital_market_gap,
        strict_portfolio_bond_gap=corrected.inherited_bond_market_gap,
        strict_portfolio_consumption_gap=strict_home_budget_gap,
        strict_portfolio_scaled_market_residual=strict_ge_residual,
        strict_lower_trade_feasibility=strict_lower,
        strict_upper_trade_feasibility=strict_upper,
        direct_lower_trade_feasibility=direct_lower,
        direct_strict_lower_trade_feasibility=direct_strict_lower,
        direct_upper_trade_feasibility=direct_upper,
        direct_strict_upper_trade_feasibility=direct_strict_upper,
        C=C,E=E,D=D,trade_price_peak=roots.peak,
        lower_trade_root=lower,upper_trade_root=upper,
        observed_feasibility=observed,
        observed_aggregate_scalar_root_gap=root_fn(w.terms[1]),
    )
end

function stage5_audit_main(root_dir::AbstractString=
    normpath(joinpath(@__DIR__,"..","results","economic_reproduction")))

    steady=solve_three_steady_states(default_parameters();verbose=false)
    boom=stage5_existing_boom(steady,root_dir)
    cases=(
        stage5_first_period_case(steady,boom,root_dir,"flexible_unleveraged";
                                 leverage_hat=0.0,fixed=false),
        stage5_first_period_case(steady,boom,root_dir,"flexible_leveraged";
                                 leverage_hat=1/16,fixed=false),
        stage5_first_period_case(steady,boom,root_dir,"fixed_leveraged";
                                 leverage_hat=1/16,fixed=true),
    )
    dest=joinpath(root_dir,"STAGE5_INITIAL_FX_ACCOUNTING_AUDIT.md")
    open(dest,"w") do io
        println(io,"# Stage 5 — Initial FX and portfolio accounting (provisional)")
        println(io)
        println(io,"This calculation loads the archived Euler-polished ",
                "Julia market roots, reconstructs household distribution ",
                "and checks **first-period aggregate portfolio accounting**.")
        println(io)
        println(io,"A positive gross-capital gap means that household ",
                "gross capital exceeds its aggregate value q_pre K_0.")
        println(io)
        for row in cases
            println(io,"## ",row.economy)
            println(io)
            println(io,"| Diagnostic | Result |")
            println(io,"|---|---:|")
            for key in (
                :maximum_scaled_residual,:terms_pre,:terms_observed,
                :observed_nominal_depreciation,
                :qK,:actual_gross_capital,:gross_capital_gap,
                :B_initial,:actual_foreign_bonds,:foreign_bond_gap,
                :pre_distribution_wealth,:pre_aggregate_target_wealth,
                :revalued_wealth,:predicted_wealth,
                :predicted_consumption_gap,:actual_consumption_gap,
                :strictly_cleared_capital,:strictly_cleared_bonds,
                :strict_portfolio_scale,:strict_portfolio_capital_gap,
                :strict_portfolio_bond_gap,:strict_portfolio_consumption_gap,
                :strict_portfolio_scaled_market_residual,
                :C,:E,:D,:trade_price_peak,
                :observed_aggregate_scalar_root_gap)
                println(io,"| ",key," | ",getproperty(row,key)," |")
            end
            println(io)
            println(io,"| Scalar trade candidate | x | Minimum disposable resources | Positive-mass affordable? |")
            println(io,"|---|---:|---:|---|")
            for (kind,branch) in (("lower",row.lower_trade_root),
                                  ("upper",row.upper_trade_root),
                                  ("observed",row.observed_feasibility))
                if branch===nothing
                    continue
                end
                println(io,"| ",kind," | ",branch.root," | ",
                        branch.minimum_resources," | ",branch.feasible," |")
                if !isempty(branch.error)
                    println(io,"  - ",kind," rejected: ",branch.error)
                end
            end
            println(io)
            println(io,"Direct (pre-interpolation) lower-root feasibility: ",
                    row.direct_lower_trade_feasibility,
                    "; strict portfolio: ",
                    row.direct_strict_lower_trade_feasibility)
            println(io)
            println(io,"Under aggregate-conserving portfolio mapping, ",
                    "lower-branch household affordability: ",
                    row.strict_lower_trade_feasibility.feasible,
                    "; upper-branch affordability: ",
                    row.strict_upper_trade_feasibility===nothing ? "no upper root" :
                        row.strict_upper_trade_feasibility.feasible,".")
            println(io)
        end
        println(io,"## Interpretation")
        println(io)
        println(io,"The scalar trade roots are **conditional on the saved aggregate K/L path**. ",
                "Multiple scalar roots do not establish multiple full 599-dimensional GE roots. ",
                "A branch may be financially infeasible for positive-mass households. ",
                "Portfolio gross capital and foreign bond totals must sum to q_pre K_0 and ",
                "B_pre respectively for the economic aggregate budget identity to hold.")
    end
    for c in cases
        println("STAGE5 ",c.economy,": capital_gap=",c.gross_capital_gap,
                ", bond_gap=",c.foreign_bond_gap,
                ", household_budget_aggregation_gap=",c.actual_consumption_gap,
                ", scalar_roots=(",
                c.lower_trade_root.root,",",
                c.upper_trade_root===nothing ? "none" : c.upper_trade_root.root,
                "), feasible_lower=",c.lower_trade_root.feasible,
                ", feasible_upper=",c.upper_trade_root===nothing ? "none" : c.upper_trade_root.feasible,
                ", strict_capital_gap=",c.strict_portfolio_capital_gap,
                ", strict_budget_gap=",c.strict_portfolio_consumption_gap,
                ", strict_599_residual_at_old_prices=",
                c.strict_portfolio_scaled_market_residual,
                ", strict_lower_feasible=",
                c.strict_lower_trade_feasibility.feasible,
                ", strict_lower_minresource=",
                c.strict_lower_trade_feasibility.minresource,
                ", strict_lower_error=",
                c.strict_lower_trade_feasibility.error,
                ", directly_revalued_original_states_lower_min=",
                c.direct_lower_trade_feasibility.minimum_resources,
                ", directly_revalued_original_states_strict_lower_min=",
                c.direct_strict_lower_trade_feasibility.minimum_resources,
                ", direct_lower_violating_mass=",
                c.direct_lower_trade_feasibility.violating_original_mass,
                ", direct_strict_lower_violating_mass=",
                c.direct_strict_lower_trade_feasibility.violating_original_mass,
                ", legacy_lower_minresource=",
                c.lower_trade_root.minimum_resources,
                ", legacy_lower_error=",c.lower_trade_root.error)
    end
    println("Saved first-period independent accounting: ",dest)
    return cases
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    isempty(ARGS) ||
        error("Usage: julia --project=. scripts/stage5_initial_fx_audit.jl")
    stage5_audit_main()
end
