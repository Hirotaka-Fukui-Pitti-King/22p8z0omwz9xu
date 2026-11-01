
"""
Stage 3 reuses the Stage 2 399-variable economic closure, but introduces
the source author's foreign-debt revaluation INSIDE EVERY candidate-path
evaluation. Its household policies and distributions are recomputed each time.
"""
function stage3_positive_trade_root(w::TransitionWorkspace{T},
                                   ctx::Stage3Context{T}) where {T}
    inp=ctx.economy
    p=inp.params
    # First-period household balance sheet after FX revaluation:
    # (q_initial*K_initial + B_initial/(x/terms_pre))*(1+r[1]).
    # The theta=1 budget and export-demand condition are otherwise
    # the identical nonlinear scalar equation in solve_transition_bis.m.
    expenditure=w.Y[1]-w.I[1]-w.adjustment[1]
    outwealth=inp.B[1]+w.q[1]*w.K[2]
    other=w.wage[1]*w.L[1]+w.profits[1]
    gross=one(T)+w.rate[1]
    function scalar_error(x::T)
        x>0 || throw(DomainError(x,"terms of trade must be positive"))
        initial_wealth=inp.q_initial*w.K[1]+inp.B_initial*ctx.terms_pre/x
        ch=p.chi*(gross*initial_wealth-outwealth+other)
        return expenditure-inp.reference.export_shift*x^(-p.theta_star)-ch
    end
    # A wide positive bracket is a diagnostic; we do not force a sign
    # change by clipping a household or its debt.
    lower=T(1e-3)*ctx.terms_pre
    upper=T(1e3)*ctx.terms_pre
    left=scalar_error(lower)
    right=scalar_error(upper)
    isfinite(left) && isfinite(right) &&
        signbit(left)!=signbit(right) ||
        throw(ArgumentError("foreign-denominated trade-price root not bracketed"))
    root=bracketed_root(scalar_error,lower,upper;tol=T(1e-11),maxiter=110)
    return root.root
end

function stage3_flexible_prices!(w::TransitionWorkspace{T},
                                 ctx::Stage3Context{T},
                                 x::AbstractVector{T}) where {T}
    inp=ctx.economy
    transition_prices!(w,inp,x)
    if ctx.foreign && inp.B_initial!=zero(T)
        price=stage3_positive_trade_root(w,ctx)
        w.terms[1]=price
        beginwealth=inp.q_initial*w.K[1]+
                    inp.B_initial*ctx.terms_pre/price
        endwealth=inp.B[1]+w.q[1]*w.K[2]
        w.home_budget[1]=inp.params.chi*
            (beginwealth*(one(T)+w.rate[1])-endwealth+
             w.wage[1]*w.L[1]+w.profits[1])
        w.home_supply[1]=w.home_budget[1]
    end
    return w
end

"""
The full grid is necessary because exchange-rate revaluation may move
positive mass below the STATIONARY borrowing-limit node. The original
backsolve_egm.m explicitly uses Agrid_fine for this. Negative-wealth
far-tail grid points with infeasible consumption are left as NaN and
must NEVER have positive distribution mass. No consumption is floored.
"""
function stage3_household_backward!(w::TransitionWorkspace{T},
                                    ctx::Stage3Context{T}) where {T}
    inp=ctx.economy
    p=inp.params
    coarse=inp.grids.coarse
    fine=inp.grids.fine
    chain=inp.chain
    nc=length(coarse); nf=length(fine); ns=length(chain.levels)
    n=length(inp.B)
    copyto!(view(w.coarse_consumption,:,:,n),
            inp.terminal.household.consumption_coarse)
    # Terminal policy uses the supplied final SS active grid; the old
    # negative tail is never needed for backward EGM (coarse only).
    fill!(view(w.fine_consumption,:,:,n),T(NaN))
    fill!(view(w.fine_savings,:,:,n),T(NaN))
    active=inp.terminal.household.consumption_fine
    terminal_grid=[v for v in fine if v>=inp.grids.borrowing_limit-T(1e-10)]
    length(terminal_grid)==size(active,1) ||
        throw(DimensionMismatch("terminal fine grid/SS policy"))
    for (i,a) in enumerate(fine)
        if a>=inp.grids.borrowing_limit-T(1e-10)
            k=searchsortedfirst(terminal_grid,a)
            k<=size(active,1) && terminal_grid[k]==a ||
                throw(ArgumentError("terminal fine nodes mismatch"))
            for j in 1:ns
                w.fine_consumption[i,j,n]=active[k,j]
                w.fine_savings[i,j,n]=inp.terminal.household.savings_fine[k,j]
            end
        end
    end
    @inbounds for t in (n-1):-1:1
        future=view(w.coarse_consumption,:,:,t+1)
        now=view(w.coarse_consumption,:,:,t)
        gross_old=one(T)+w.rate[t]
        gross_new=one(T)+w.rate[t+1]
        m=consumption_cost(w.terms[t],p)
        for j in 1:ns
            income=w.wage[t]*w.L[t]*chain.levels[j]+w.profits[t]
            for i in 1:nc
                emu=zero(T)
                for next_s in 1:ns
                    emu+=chain.transition[j,next_s]*future[i,next_s]^(-p.gamma)
                end
                cc=(inp.initial.beta*gross_new*emu)^(-inv(p.gamma))
                w.endogenous_consumption[i]=cc
                w.endogenous_asset[i]=(m*cc+coarse[i]-income)/gross_old
            end
            all(diff(w.endogenous_asset) .> 0) ||
                throw(ArgumentError("non-monotone endogenous wealth at time=$t, s=$j"))
            pchip_slopes!(w.interpolation_derivative,w.endogenous_asset,
                          w.endogenous_consumption)
            minendo=w.endogenous_asset[1]
            for i in 1:nc
                a=coarse[i]
                c=a<=minendo ?
                  (gross_old*a+income-inp.grids.borrowing_limit)/m :
                  pchip_evaluate(w.endogenous_asset,w.endogenous_consumption,
                                 w.interpolation_derivative,a)
                c>0 && isfinite(c) ||
                    throw(DomainError(c,"invalid feasible coarse EGM consumption"))
                now[i,j]=c
            end
            for i in 1:nf
                a=fine[i]
                if a<=minendo
                    saving=inp.grids.borrowing_limit
                    cons=(gross_old*a+income-saving)/m
                else
                    cons=pchip_evaluate(w.endogenous_asset,
                                       w.endogenous_consumption,
                                       w.interpolation_derivative,a)
                    saving=gross_old*a+income-m*cons
                end
                if !(isfinite(cons) && cons>0 &&
                     isfinite(saving) && saving>=inp.grids.borrowing_limit)
                    # Negative pre-constraint nodes not carrying positive
                    # probability are unavailable; validate after
                    # revaluation/forward transport, never silently clip.
                    if a>=inp.grids.borrowing_limit
                        throw(DomainError((a,cons,saving),
                                          "infeasible nonnegative-domain EGM"))
                    end
                    cons=T(NaN)
                    saving=T(NaN)
                end
                w.fine_consumption[i,j,t]=cons
                w.fine_savings[i,j,t]=saving
            end
        end
    end
    return w
end

function stage3_distribution_forward!(w::TransitionWorkspace{T},
                                      ctx::Stage3Context{T}) where {T}
    inp=ctx.economy
    grid=inp.grids.fine
    init=view(w.distribution,:,:,1)
    revalue_foreign_debt!(init,ctx.mu_pre,grid,w.K[1],
                          inp.q_initial,inp.B_initial,ctx.terms_pre,
                          w.terms[1],ctx.leverage_hat;foreign=ctx.foreign)
    ns=length(inp.chain.levels); nf=length(grid); n=length(inp.B)
    @inbounds for t in 1:n
        mu=view(w.distribution,:,:,t)
        s=zero(T); ch=zero(T)
        for j in 1:ns,i in 1:nf
            mass=mu[i,j]
            mass==0 && continue
            mass>0 || throw(DomainError(mass,"negative household mass"))
            c=w.fine_consumption[i,j,t]
            a=w.fine_savings[i,j,t]
            isfinite(c) && c>0 && isfinite(a) &&
                a>=inp.grids.borrowing_limit ||
                throw(DomainError((t,grid[i],j,c,a),
                                   "positive-mass household infeasible after FX shock"))
            ch+=mass*c
            s+=mass*a
        end
        w.aggregate_home[t]=ch
        w.aggregate_savings[t]=s
        w.distribution_error[t]=abs(sum(mu)-one(T))
        w.distribution_error[t]<T(1e-8) ||
            throw(DomainError(w.distribution_error[t],"distribution mass loss"))
        if t<n
            distribution_step!(view(w.distribution,:,:,t+1),
                               mu,view(w.fine_savings,:,:,t),
                               grid,inp.chain.transition)
        end
    end
    return w
end

function stage3_market_residual!(res::AbstractVector{T},
                                 w::TransitionWorkspace{T},
                                 ctx::Stage3Context{T};
                                 fixed::Bool=false) where {T}
    inp=ctx.economy
    p=inp.params; n=length(inp.B)
    expected=fixed ? 3n-1 : 2n-1
    length(res)==expected || throw(DimensionMismatch("Stage3 GE equation count"))
    scaleK=max(one(T),inp.initial.K)
    scaleL=max(one(T),inp.initial.L)
    @inbounds for t in 1:n-1
        res[t]=((w.aggregate_savings[t]-inp.B[t])/w.q[t]-w.K[t+1])/scaleK
    end
    unionfactor=p.chi*(p.epsilon_w-1)/p.epsilon_w
    @inbounds for t in 1:n
        c=t<n ? w.aggregate_home[t] : w.aggregate_home[t-1]
        c>0 || throw(DomainError(c,"negative aggregate home consumption"))
        w.labor_implied[t]=unionfactor*w.wage[t]/
            (inp.reference.chi_dis*c)
        res[n-1+t]=(w.labor_implied[t]-w.L[t])/scaleL
        if fixed
            inflation=t==1 ? w.terms[1]/ctx.terms_pre :
                              w.terms[t]/w.terms[t-1]
            nominal_cost=ctx.theta_nr/T(2)*(inflation-one(T))^2
            if t<n
                goods=w.Y[t]-w.I[t]-w.adjustment[t]-
                    nominal_cost-inp.reference.export_shift*
                    w.terms[t]^(-p.theta_star)
            else
                # Replicate final-period investment/boundary timing.
                goods=w.Y[t]-w.I[t-1]-nominal_cost-
                    inp.reference.export_shift*w.terms[t]^(-p.theta_star)-
                    ctx.economy.phi_ac/T(2)*w.K[t]*
                    ((w.I[t-1]-p.delta*w.K[t])/w.K[t])^2
            end
            res[2n-1+t]=(goods-c)/max(one(T),inp.reference.Y)
        end
    end
    all(isfinite,res) || throw(ArgumentError("nonfinite Stage 3 GE residual"))
    return res
end

function stage3_residual!(res::AbstractVector{T}, x::AbstractVector{T},
                          ctx::Stage3Context{T},
                          w::TransitionWorkspace{T}) where {T}
    w.function_evals+=1
    stage3_flexible_prices!(w,ctx,x)
    stage3_household_backward!(w,ctx)
    stage3_distribution_forward!(w,ctx)
    return stage3_market_residual!(res,w,ctx)
end

"""
Dense damped finite-difference Broyden/Newton solves the exact heterogeneous
household general equilibrium each trial. Shared between the 399-variable
flexible and 599-variable fixed models; the economic residuals differ.
"""
function stage3_nonlinear_solve(ctx::Stage3Context{T},x0::Vector{T};
                                fixed::Bool=false,tolerance::T=T(1e-6),
                                maxiter::Int=70,refresh::Int=15,
                                verbose::Bool=true) where {T}
    n=length(x0)
    expected=fixed ? 3*length(ctx.economy.B)-1 :
                     2*length(ctx.economy.B)-1
    n==expected || throw(DimensionMismatch("Stage3 unknown dimensions"))
    work=TransitionWorkspace(ctx.economy)
    evaluation! = fixed ?
        ((res,x)->stage3_fixed_residual!(res,x,ctx,work)) :
        ((res,x)->stage3_residual!(res,x,ctx,work))
    x=copy(x0); f=zeros(T,n); temp=similar(f); trial=similar(x)
    J=Matrix{T}(undef,n,n)
    evaluation!(f,x)
    evaluations=1
    for iter in 1:maxiter
        error_max=maximum(abs,f)
        verbose && println("STAGE3 ",fixed ? "FIXED" : "FLEX",
                            " iteration=",iter-1," scaled=",error_max,
                            " evaluations=",work.function_evals)
        if error_max<=tolerance
            evaluation!(f,x)
            return TransitionSolveResult(ctx.economy,work,x,copy(f),
                        iter-1,work.function_evals,maximum(abs,f),true)
        end
        if iter==1 || (iter-1)%refresh==0
            @inbounds for j in 1:n
                copyto!(trial,x)
                dx=T(1e-5)*max(one(T),abs(x[j]))
                trial[j]+=dx
                evaluation!(temp,trial)
                for i in 1:n
                    J[i,j]=(temp[i]-f[i])/dx
                end
            end
            evaluation!(f,x)
        end
        direction=-(J\f)
        all(isfinite,direction) || error("nonfinite Stage3 Newton direction")
        accepted=false
        alpha=one(T)
        for _ in 1:18
            @inbounds for j in 1:n
                trial[j]=x[j]+alpha*direction[j]
            end
            valid=true
            try
                evaluation!(temp,trial)
            catch err
                if err isa DomainError || err isa ArgumentError
                    valid=false
                else
                    rethrow()
                end
            end
            if valid && norm(temp)<(one(T)-T(1e-4)*alpha)*norm(f)
                accepted=true
                break
            end
            alpha/=T(2)
        end
        accepted || error("Stage3 line search failed at iteration=$iter; max=$(maximum(abs,f))")
        step=alpha.*direction
        correction=(temp.-f).-J*step
        denom=dot(step,step)
        denom>0 || error("vanishing Stage3 nonlinear step")
        @inbounds for j in 1:n, i in 1:n
            J[i,j]+=correction[i]*step[j]/denom
        end
        copyto!(x,trial)
        copyto!(f,temp)
    end
    error("Stage3 GE not converged: max=$(maximum(abs,f)), evaluations=$(work.function_evals)")
end

function solve_stage3_flexible(ctx::Stage3Context{T};
                               x0::AbstractVector{T}=initial_transition_guess(ctx.economy),
                               kwargs...) where {T}
    return stage3_nonlinear_solve(ctx,Vector{T}(x0);fixed=false,kwargs...)
end
