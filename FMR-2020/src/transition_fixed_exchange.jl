
"""
Fixed nominal exchange rate, foreign-denominated balance sheets:
the original solve_transition_fixed.m is a 3T-1 unknown problem.
K[2:T], L[1:T], and terms[1:T] are endogenous. The last T residuals are (T-1) actual domestic goods-market equations
and the stationary continuation terms-of-trade boundary. Firm real marginal costs are solved
backwards from the original nonlinear Rotemberg price-setting equation.
"""
function stage3_fixed_prices!(w::TransitionWorkspace{T},
                              ctx::Stage3Context{T},
                              x::AbstractVector{T}) where {T}
    inp=ctx.economy
    p=inp.params
    n=length(inp.B)
    length(x)==3*n-1 ||
        throw(DimensionMismatch("fixed-rate transition needs 3T-1 unknowns"))
    w.K[1]=inp.initial.K
    @inbounds for t in 2:n
        w.K[t]=x[t-1]
    end
    @inbounds for t in 1:n
        w.L[t]=x[n-1+t]
        w.terms[t]=x[2*n-1+t]
        w.K[t]>0 && w.L[t]>0 && w.terms[t]>0 ||
            throw(DomainError((w.K[t],w.L[t],w.terms[t]),
                               "positive capital, labor and terms required"))
        w.Y[t]=inp.z[t]*w.K[t]^p.alpha*w.L[t]^(one(T)-p.alpha)
        w.I[t]=t<n ? w.K[t+1]-(one(T)-p.delta)*w.K[t] :
                        p.delta*w.K[t]
        if t<n
            gap=w.I[t]/w.K[t]-p.delta
            w.q[t]=one(T)+inp.phi_ac*gap
            w.adjustment[t]=inp.phi_ac/T(2)*w.K[t]*gap^2
        else
            w.q[t]=one(T)
            w.adjustment[t]=zero(T)
        end
        w.q[t]>0 || throw(DomainError(w.q[t],"positive capital shadow price"))
    end
    gross=Vector{T}(undef,n)
    infl=Vector{T}(undef,n)
    @inbounds for t in 1:n
        infl[t]=t==1 ? w.terms[t]/ctx.terms_pre :
                        w.terms[t]/w.terms[t-1]
        infl[t]>0 || throw(DomainError(infl[t],"positive gross inflation"))
    end
    markup=(p.epsilon_f-one(T))/p.epsilon_f
    mc=Vector{T}(undef,n)
    mc[n]=markup+ctx.theta_nr/p.epsilon_f*
          (infl[n]-one(T))*infl[n]
    @inbounds for t in n:-1:1
        if t<n
            mc[t]=markup+ctx.theta_nr/p.epsilon_f*
                    (infl[t]-one(T))*infl[t]-
                  ctx.theta_nr/p.epsilon_f*
                    (infl[t+1]-one(T))*infl[t+1]*
                    w.Y[t+1]/(gross[t+1]*w.Y[t])
        end
        mc[t]>0 || throw(DomainError(mc[t],
                                    "negative real marginal cost at period $t"))
        w.wage[t]=mc[t]*(one(T)-p.alpha)*w.Y[t]/w.L[t]
        w.rental[t]=mc[t]*p.alpha*w.Y[t]/w.K[t]
        qprev=t==1 ? inp.q_initial : w.q[t-1]
        w.rate[t]=(w.rental[t]+(one(T)-p.delta)*w.q[t])/qprev-one(T)
        gross[t]=one(T)+w.rate[t]
        gross[t]>0 || throw(DomainError(w.rate[t],
                                       "nonpositive gross household return"))
        # Original firms' profits: Y - rK*K - wage*L
        # - theta_nr/2*(gross inflation-1)^2.
        nominal_cost=ctx.theta_nr/T(2)*(infl[t]-one(T))^2
        w.profits[t]=w.Y[t]-w.rental[t]*w.K[t]-
                      w.wage[t]*w.L[t]-nominal_cost
        if t<n
            w.profits[t]+=(w.q[t]-one(T))*w.I[t]-w.adjustment[t]
        end
    end
    @inbounds for t in 1:n
        beginning=t==1 ?
            inp.q_initial*w.K[1]+
               (ctx.foreign ? inp.B_initial*ctx.terms_pre/w.terms[1] :
                              inp.B_initial) :
            inp.B[t-1]+w.q[t-1]*w.K[t]
        ending=t<n ? inp.B[t]+w.q[t]*w.K[t+1] :
                     inp.B[t]+inp.terminal.K
        w.home_budget[t]=p.chi*
            (beginning*(one(T)+w.rate[t])-ending+
             w.wage[t]*w.L[t]+w.profits[t])
        w.home_supply[t]=w.home_budget[t]
    end
    return w
end

function stage3_fixed_residual!(res::AbstractVector{T},
                                x::AbstractVector{T},
                                ctx::Stage3Context{T},
                                w::TransitionWorkspace{T};
                                polish_households::Bool=false,
                                strict_portfolio::Bool=false) where {T}
    w.function_evals+=1
    stage3_fixed_prices!(w,ctx,x)
    stage3_household_backward!(w,ctx)
    if polish_households
        stage3_euler_polish_households!(w,ctx)
    end
    stage3_distribution_forward!(w,ctx;strict_portfolio=strict_portfolio)
    return stage3_market_residual!(res,w,ctx;fixed=true)
end

function solve_stage3_fixed(ctx::Stage3Context{T},
                            prior::TransitionSolveResult{T};
                            kwargs...) where {T}
    n=length(ctx.economy.B)
    prior.converged && length(prior.x)==2*n-1 ||
        throw(ArgumentError("requires solved flexible sudden-stop path"))
    x=vcat(prior.x,prior.workspace.terms)
    return stage3_nonlinear_solve(ctx,x;fixed=true,kwargs...)
end

"""
Independent Stage 3 checks: the expanded inherited negative wealth
support is allowed only if household consumption and saving are valid
where the distribution has POSITIVE probability. Fixed-regime goods clearing uses actual cross-sectional C_H at all
T-1 economic transition dates; date T is a stationary continuation boundary.
"""
function stage3_diagnostics(sol::TransitionSolveResult{T},
                            ctx::Stage3Context{T};fixed::Bool=false,
                            explicit_flexible::Bool=false) where {T}
    sol.converged || throw(ArgumentError("Stage3 unaccepted GE candidate"))
    w=sol.workspace; inp=ctx.economy
    n=length(inp.B)
    size(w.distribution,3)==n ||
        throw(DimensionMismatch("distribution horizon"))
    maxasset=zero(T); maxlabor=zero(T); maxgoods=zero(T)
    maxbudget=zero(T); maxterminalbudget=zero(T); maxmass=zero(T)
    maxphillips=zero(T)
    maxeuler=zero(T); weul=zero(T); wmass=zero(T)
    maxeuler_first40=zero(T); maxeuler_last20=zero(T)
    weighted_terminal_budget=zero(T); terminal_budget_mass=zero(T)
    euler_mass_above_1e4=zero(T)
    euler_states_above_1e4=0
    maxeuler_state=(period=0,asset_index=0,productivity_index=0,
                    mass=zero(T),asset=zero(T),savings=zero(T),
                    binding=false,relative_gap=zero(T))
    grid=inp.grids.fine
    slopes=Matrix{T}(undef,length(inp.grids.coarse),length(inp.chain.levels))
    @inbounds for t in 1:n
        if t<n
            future=view(w.coarse_consumption,:,:,t+1)
            for jp in axes(future,2)
                pchip_slopes!(view(slopes,:,jp),inp.grids.coarse,
                              view(future,:,jp))
            end
        end
        mu=view(w.distribution,:,:,t)
        maxmass=max(maxmass,abs(sum(mu)-one(T)))
        if t<n
            maxlabor=max(maxlabor,abs(w.labor_implied[t]-w.L[t]))
        end
        if t<n
            maxasset=max(maxasset,
                 abs((w.aggregate_savings[t]-inp.B[t])/w.q[t]-w.K[t+1]))
        end
        infl=t==1 ? w.terms[t]/ctx.terms_pre :
                    w.terms[t]/w.terms[t-1]
        nominalcost=fixed ? ctx.theta_nr/T(2)*(infl-one(T))^2 : zero(T)
        if t<n
            supply=w.Y[t]-w.I[t]-w.adjustment[t]-nominalcost-
                    inp.reference.export_shift*w.terms[t]^(-inp.params.theta_star)
            maxgoods=max(maxgoods,abs(supply-w.aggregate_home[t]))
        end
        for j in axes(mu,2),i in axes(mu,1)
            mass=mu[i,j]
            mass==0 && continue
            a=grid[i]
            c=w.fine_consumption[i,j,t]
            saving=w.fine_savings[i,j,t]
            mass>0 && isfinite(c) && c>0 && isfinite(saving) ||
                throw(DomainError((mass,c,saving),"invalid occupied household"))
            if t<n
                m=consumption_cost(w.terms[t],inp.params)
                income=w.wage[t]*w.L[t]*inp.chain.levels[j]+w.profits[t]
                current_gap=(one(T)+w.rate[t])*a+income-
                            m*c-saving
                maxbudget=max(maxbudget,abs(current_gap))
                # Independent Euler evaluation on coarse next-date policy
                future=view(w.coarse_consumption,:,:,t+1)
                emu=zero(T)
                for jp in axes(future,2)
                    nextc=pchip_evaluate(inp.grids.coarse,
                                         view(future,:,jp),
                                         view(slopes,:,jp),saving)
                    nextc>0 || throw(DomainError(nextc,"Euler future consumption"))
                    emu+=inp.chain.transition[j,jp]*nextc^(-inp.params.gamma)
                end
                lhs=c^(-inp.params.gamma)
                rhs=inp.initial.beta*(one(T)+w.rate[t+1])*emu
                gap=(lhs-rhs)/lhs
                violation=saving<=inp.grids.borrowing_limit+T(1e-8) ?
                         max(zero(T),-gap) : abs(gap)
                if violation>maxeuler
                    maxeuler=violation
                    maxeuler_state=(period=t,asset_index=i,
                                    productivity_index=j,mass=mass,
                                    asset=a,savings=saving,
                                    binding=saving<=inp.grids.borrowing_limit+T(1e-8),
                                    relative_gap=gap)
                end
                if t<=min(n-1,40)
                    maxeuler_first40=max(maxeuler_first40,violation)
                end
                if t>=max(1,n-20)
                    maxeuler_last20=max(maxeuler_last20,violation)
                end
                if violation>T(1e-4)
                    euler_mass_above_1e4+=mass
                    euler_states_above_1e4+=1
                end
                weul+=mass*violation
                wmass+=mass
            else
                # Terminal c and a' come from the stationary continuation
                # problem. Their budget feasibility at actual period-T
                # prices is an independent horizon-truncation check,
                # not a market-clearing equation to be hidden by using C[T-1].
                m=consumption_cost(w.terms[t],inp.params)
                income=w.wage[t]*w.L[t]*inp.chain.levels[j]+w.profits[t]
                budget=(one(T)+w.rate[t])*a+income-m*c-saving
                maxterminalbudget=max(maxterminalbudget,abs(budget))
                weighted_terminal_budget+=mass*abs(budget)
                terminal_budget_mass+=mass
            end
        end
    end
    # A direct one-dimensional constrained household optimal-choice check
    # at the worst occupied Euler state. This compares the EGM interpolated
    # savings policy with the exact numerical best response to the SAME
    # continuous continuation policy, without a source MATLAB reference.
    worst_best_response=nothing
    if maxeuler_state.period>0
        t=maxeuler_state.period
        j=maxeuler_state.productivity_index
        i=maxeuler_state.asset_index
        future=view(w.coarse_consumption,:,:,t+1)
        income=w.wage[t]*w.L[t]*inp.chain.levels[j]+w.profits[t]
        gross=one(T)+w.rate[t]
        resources=gross*grid[i]+income
        candidate=household_best_response_saving(
            inp.grids.coarse,future,view(inp.chain.transition,j,:),
            inp.params.gamma,inp.initial.beta,one(T)+w.rate[t+1],
            resources,consumption_cost(w.terms[t],inp.params),
            inp.grids.borrowing_limit)
        worst_best_response=(
            numerical_best_savings=candidate.savings,
            egm_savings=w.fine_savings[i,j,t],
            absolute_savings_difference=abs(candidate.savings-w.fine_savings[i,j,t]),
            numerical_best_consumption=candidate.consumption,
            egm_consumption=w.fine_consumption[i,j,t],
            binding=candidate.binding,
            direct_euler_surplus=candidate.euler_surplus,
            bisection_iterations=candidate.iterations)
    end
    if fixed
        # Reconstruct every nonlinear NKPC equation independently from
        # realized real MC=(wage*L)/[(1-alpha)Y].
        for t in 1:n
            pi=t==1 ? w.terms[t]/ctx.terms_pre :
                      w.terms[t]/w.terms[t-1]
            markup=(inp.params.epsilon_f-one(T))/inp.params.epsilon_f
            desired=markup+ctx.theta_nr/inp.params.epsilon_f*(pi-one(T))*pi
            if t<n
                nextpi=w.terms[t+1]/w.terms[t]
                desired-=ctx.theta_nr/inp.params.epsilon_f*
                         (nextpi-one(T))*nextpi*
                         w.Y[t+1]/((one(T)+w.rate[t+1])*w.Y[t])
            end
            actual=w.wage[t]*w.L[t]/
                   ((one(T)-inp.params.alpha)*w.Y[t])
            maxphillips=max(maxphillips,abs(actual-desired))
        end
    end
    return (max_asset_market_level=maxasset,
            max_labor_market_level=maxlabor,
            max_goods_market_level=maxgoods,
            max_household_budget=maxbudget,
            max_terminal_household_budget=maxterminalbudget,
            weighted_terminal_household_budget=weighted_terminal_budget/
                max(terminal_budget_mass,eps(T)),
            max_mass_deviation=maxmass,
            max_nonlinear_phillips=maxphillips,
            max_dynamic_euler=maxeuler,
            max_dynamic_euler_first40=maxeuler_first40,
            max_dynamic_euler_last20=maxeuler_last20,
            max_dynamic_euler_state=maxeuler_state,
            worst_state_direct_best_response=worst_best_response,
            euler_mass_above_1e4=euler_mass_above_1e4/wmass,
            euler_states_above_1e4=euler_states_above_1e4,
            weighted_mean_dynamic_euler=weul/wmass,
            maximum_scaled_solver_residual=sol.maximum_scaled_residual,
            iterations=sol.iterations,evaluations=sol.evaluations,
            terminal_K_gap=abs(w.K[end]-inp.terminal.K),
            terminal_L_gap=abs(w.L[end]-inp.terminal.L),
            terminal_terms_gap=abs(w.terms[end]-inp.terminal.terms),
            terminal_interest_gap=abs(w.rate[end]-inp.terminal.r),
            terminal_wage_gap=abs(w.wage[end]-inp.terminal.wage),
            terminal_profit_gap=abs(w.profits[end]-inp.terminal.profits))
end


"""
    stage3_euler_polish_households!(w, ctx; trigger=1e-5)

Correct the *fixed-price* Stage3 household policy by solving the
one-period Euler/KKT equation directly only at coarse/fine asset
states where inverse-EGM interpolation has a material discrepancy.
Backwards through time: corrected coarse policies become next-period
continuation values for earlier decisions. Keeps the original stationary
terminal boundary and all prices/aggregate guesses unchanged.

DO NOT treat the resulting policies as a market-clearing equilibrium.
After polishing, the caller MUST propagate distribution forward and
recompute the full aggregate residual at the same K/L/terms candidate.
If that residual changes, the aggregate GE must be resolved using the
polished household map. Intended first as a one-pass numerical audit.
"""
function stage3_euler_polish_households!(
    w::TransitionWorkspace{T},ctx::Stage3Context{T};
    trigger::T=T(1e-5),target::T=T(1e-9)) where {T<:AbstractFloat}

    inp=ctx.economy
    coarse=inp.grids.coarse
    fine=inp.grids.fine
    chain=inp.chain
    p=inp.params
    n=length(inp.B)
    nc=length(coarse)
    ns=length(chain.levels)
    slopes=Matrix{T}(undef,nc,ns)
    improved_coarse=0
    improved_fine=0
    offgrid_fine=0
    maximum_savings_update=zero(T)
    max_post_gap=zero(T)
    worst_post=(period=0,asset_index=0,productivity_index=0)
    @inbounds for t in (n-1):-1:1
        future=view(w.coarse_consumption,:,:,t+1)
        now=view(w.coarse_consumption,:,:,t)
        cnow=view(w.fine_consumption,:,:,t)
        snow=view(w.fine_savings,:,:,t)
        for next_j in 1:ns
            pchip_slopes!(view(slopes,:,next_j),coarse,
                          view(future,:,next_j))
        end
        R=one(T)+w.rate[t]
        Rnext=one(T)+w.rate[t+1]
        m=consumption_cost(w.terms[t],p)
        for j in 1:ns
            probs=view(chain.transition,j,:)
            income=w.wage[t]*w.L[t]*chain.levels[j]+w.profits[t]
            for i in eachindex(coarse)
                resources=R*coarse[i]+income
                sav=resources-m*now[i,j]
                # Only machine-scale budget roundoff can be normalized.
                if sav<inp.grids.borrowing_limit &&
                   sav>=inp.grids.borrowing_limit-T(1e-12)
                    sav=inp.grids.borrowing_limit
                end
                result=euler_refine_one_saving(
                    coarse,future,slopes,probs,p.gamma,inp.initial.beta,
                    Rnext,resources,m,inp.grids.borrowing_limit,sav;
                    trigger=trigger,target=target)
                if result.refined
                    improved_coarse+=1
                    now[i,j]=result.consumption
                end
            end
            for i in eachindex(fine)
                sav=snow[i,j]
                if !isfinite(sav) || !isfinite(cnow[i,j])
                    continue
                end
                resources=R*fine[i]+income
                result=euler_refine_one_saving(
                    coarse,future,slopes,probs,p.gamma,inp.initial.beta,
                    Rnext,resources,m,inp.grids.borrowing_limit,sav;
                    trigger=trigger,target=target)
                if result.offgrid
                    offgrid_fine+=1
                    continue
                end
                if result.refined
                    improved_fine+=1
                    maximum_savings_update=max(maximum_savings_update,
                                               abs(result.saving-sav))
                    snow[i,j]=result.saving
                    cnow[i,j]=result.consumption
                end
                # At a binding constraint, positive Euler surplus is
                # compatible with optimality (KKT), not a violation.
                violation=result.saving<=inp.grids.borrowing_limit+T(1e-8) ?
                    max(zero(T),-result.relative_gap) :
                    abs(result.relative_gap)
                if violation>max_post_gap
                    max_post_gap=violation
                    worst_post=(period=t,asset_index=i,
                                productivity_index=j)
                end
            end
        end
    end
    return (refined_coarse=improved_coarse,
            refined_fine=improved_fine,
            offgrid_fine_states=offgrid_fine,
            maximum_savings_update=maximum_savings_update,
            maximum_euler_violation_on_supported_grid=max_post_gap,
            worst_post_state=worst_post)
end
