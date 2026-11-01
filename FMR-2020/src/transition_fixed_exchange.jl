
"""
Fixed nominal exchange rate, foreign-denominated balance sheets:
the original solve_transition_fixed.m is a 3T-1 unknown problem.
K[2:T], L[1:T], and terms[1:T] are endogenous. The last T residuals
are domestic goods clearing. Firm real marginal costs are solved
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
                                w::TransitionWorkspace{T}) where {T}
    w.function_evals+=1
    stage3_fixed_prices!(w,ctx,x)
    stage3_household_backward!(w,ctx)
    stage3_distribution_forward!(w,ctx)
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
where the distribution has POSITIVE probability. Fixed-regime goods
clearing uses actual cross-sectional C_H in ALL T equations.
"""
function stage3_diagnostics(sol::TransitionSolveResult{T},
                            ctx::Stage3Context{T};fixed::Bool=false) where {T}
    sol.converged || throw(ArgumentError("Stage3 unaccepted GE candidate"))
    w=sol.workspace; inp=ctx.economy
    n=length(inp.B)
    size(w.distribution,3)==n ||
        throw(DimensionMismatch("distribution horizon"))
    maxasset=zero(T); maxlabor=zero(T); maxgoods=zero(T)
    maxbudget=zero(T); maxmass=zero(T)
    maxphillips=zero(T)
    maxeuler=zero(T); weul=zero(T); wmass=zero(T)
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
        maxlabor=max(maxlabor,abs(w.labor_implied[t]-w.L[t]))
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
        elseif fixed
            cost=inp.phi_ac/T(2)*w.K[t]*
                     ((w.I[t-1]-inp.params.delta*w.K[t])/w.K[t])^2
            supply=w.Y[t]-w.I[t-1]-cost-nominalcost-
                    inp.reference.export_shift*w.terms[t]^(-inp.params.theta_star)
            maxgoods=max(maxgoods,abs(supply-w.aggregate_home[t-1]))
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
                maxeuler=max(maxeuler,violation)
                weul+=mass*violation
                wmass+=mass
            end
        end
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
            max_mass_deviation=maxmass,
            max_nonlinear_phillips=maxphillips,
            max_dynamic_euler=maxeuler,
            weighted_mean_dynamic_euler=weul/wmass,
            maximum_scaled_solver_residual=sol.maximum_scaled_residual,
            iterations=sol.iterations,evaluations=sol.evaluations,
            terminal_K_gap=abs(w.K[end]-inp.terminal.K),
            terminal_L_gap=abs(w.L[end]-inp.terminal.L))
end
