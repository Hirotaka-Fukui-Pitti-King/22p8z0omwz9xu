
# Stage 3 reuses the Stage 2 399-variable economic closure, but introduces
# the source author's foreign-debt revaluation INSIDE EVERY candidate-path
# evaluation. Its household policies and distributions are recomputed each time.
#
"""
For negative opening foreign bonds, the source trade-price equation is

    f(x) = C + E/x - D/x^theta_star = 0, E > 0, D > 0.

For theta_star > 1 its only interior stationary point is a MAXIMUM at
x_peak=(theta_star*D/E)^(1/(theta_star-1)). If C < 0 it can have
TWO positive roots despite f(lower)<0 and f(upper)<0.

The original MATLAB solve_transition_bis.m calls fsolve(aggr,0.8) in
its final iteration. Preserve the LOW-price branch (the one reached
from 0.8 in the observed PC initial candidate), rather than
incorrectly testing opposite signs at arbitrary far-apart endpoints.
Do not change aggregate paths or clip households to manufacture a root.
"""
function foreign_debt_trade_root(C::T,E::T,D::T,theta_star::T) where {T<:AbstractFloat}
    isfinite(C) && isfinite(E) && isfinite(D) && isfinite(theta_star) ||
        throw(DomainError((C,E,D,theta_star),"nonfinite scalar trade coefficients"))
    E>zero(T) && D>zero(T) && theta_star>one(T) ||
        throw(DomainError((E,D,theta_star),
                          "root-branch formula requires foreign debt E>0, D>0, theta_star>1"))
    xpeak=(theta_star*D/E)^(inv(theta_star-one(T)))
    f(x)=C+E/x-D*x^(-theta_star)
    fpeak=f(xpeak)
    isfinite(xpeak) && xpeak>zero(T) && isfinite(fpeak) && fpeak>zero(T) ||
        throw(DomainError((C,E,D,xpeak,fpeak),
            "no two-sided positive lower-branch trade root at the proposed aggregate path"))
    low=xpeak/T(2)
    flow=f(low)
    for _ in 1:100
        isfinite(flow) && flow<zero(T) && break
        low/=T(2)
        flow=f(low)
    end
    isfinite(flow) && flow<zero(T) ||
        throw(DomainError((low,flow),"failed to bracket positive low-price branch"))
    root=bracketed_root(f,low,xpeak;tol=T(1e-11),maxiter=110)
    isfinite(root.root) && zero(T)<root.root<xpeak ||
        throw(DomainError(root.root,"invalid source low-price root"))
    abs(f(root.root))<=T(1e-9)*
                          max(one(T),abs(C),abs(E/root.root),abs(D/root.root^theta_star)) ||
        throw(ErrorException("foreign-debt trade root has unacceptable equation residual"))
    return root.root
end


"""
Compute all strictly positive roots of the original scalar trade equation.
The low root alone is NOT generally a feasible heterogeneous-household GE:
after FX debt revaluation, occupied states may violate the household budget.
"""
function foreign_debt_trade_roots(C::T,E::T,D::T,theta_star::T) where {T<:AbstractFloat}
    low=foreign_debt_trade_root(C,E,D,theta_star)
    peak=(theta_star*D/E)^(inv(theta_star-one(T)))
    high=nothing
    if C<zero(T)
        f(x)=C+E/x-D*x^(-theta_star)
        upper=T(2)*peak
        for _ in 1:100
            f(upper)<zero(T) && break
            upper*=T(2)
        end
        isfinite(upper) && f(upper)<zero(T) ||
            throw(DomainError(upper,"unable to bracket upper positive trade root"))
        root=bracketed_root(f,peak,upper;tol=T(1e-11),maxiter=110)
        high=root.root
        high>peak || throw(DomainError(high,"upper trade root not above peak"))
    end
    return (low=low,high=high,peak=peak)
end

"""
Necessary feasibility test, before backward EGM, on the exact positively
weighted grid states created by the original distribution mass-splitting
scheme. It rejects the author MATLAB's numerical consumption floor when that
floor would give a household resources below its exogenous borrowing limit.
No zero-mass exception, borrowing-limit relaxation, or probability clipping.
"""
function initial_revaluation_budget_support!(w::TransitionWorkspace{T},
                                              ctx::Stage3Context{T},
                                              terms::T) where {T}
    inp=ctx.economy
    p=inp.params
    grid=inp.grids.fine
    mu=view(w.distribution,:,:,1)
    revalue_foreign_debt!(mu,ctx.mu_pre,grid,w.K[1],
                          inp.q_initial,inp.B_initial,
                          ctx.terms_pre,terms,ctx.leverage_hat;
                          foreign=ctx.foreign)
    gross=one(T)+w.rate[1]
    minresource=T(Inf)
    minstate=(0,0)
    for j in axes(mu,2),i in axes(mu,1)
        mu[i,j]>zero(T) || continue
        earnings=w.wage[1]*w.L[1]*inp.chain.levels[j]+w.profits[1]
        disposable=gross*grid[i]+earnings-inp.grids.borrowing_limit
        if disposable<minresource
            minresource=disposable
            minstate=(i,j)
        end
    end
    return (feasible=isfinite(minresource)&&minresource>zero(T),
            minresource=minresource,minstate=minstate,
            occupied_mass=sum(mu))
end

function stage3_positive_trade_root(w::TransitionWorkspace{T},
                                   ctx::Stage3Context{T}) where {T}
    inp=ctx.economy
    p=inp.params
    expenditure=w.Y[1]-w.I[1]-w.adjustment[1]
    outwealth=inp.B[1]+w.q[1]*w.K[2]
    other=w.wage[1]*w.L[1]+w.profits[1]
    gross=one(T)+w.rate[1]
    C=expenditure-p.chi*(gross*inp.q_initial*w.K[1]-outwealth+other)
    E=-p.chi*gross*inp.B_initial*ctx.terms_pre
    D=inp.reference.export_shift
    E>zero(T) ||
        throw(DomainError(E,"negative inherited foreign bonds required for this root map"))
    roots=foreign_debt_trade_roots(C,E,D,p.theta_star)
    # MATLAB fsolve(aggr,0.8) may return the lower root, but it also
    # imposes c=max(c_budget,1e-5) even on households whose budgets
    # CANNOT fund c=1e-5. The correct GE must pass household feasibility.
    # Prefer the author's lower branch only if it is truly feasible.
    # An out-of-grid positive-mass node rejects THAT root only; it
    # must not prevent evaluating the other economically possible root.
    function affordable(candidate::T)
        try
            return initial_revaluation_budget_support!(w,ctx,candidate)
        catch err
            if err isa DomainError
                return (feasible=false,minresource=-T(Inf),
                        minstate=(0,0),occupied_mass=zero(T))
            end
            rethrow()
        end
    end
    lowcheck=affordable(roots.low)
    highcheck=roots.high===nothing ? nothing :
              affordable(roots.high)
    selected=if lowcheck.feasible
        roots.low
    elseif highcheck!==nothing && highcheck.feasible
        roots.high
    else
        throw(DomainError((roots.low,roots.high,lowcheck,
                           highcheck),
            "neither trade-price branch supports positive-mass household budgets"))
    end
    if w.function_evals==1
        println("STAGE3 TRADE ROOTS low=",roots.low,
                " high=",roots.high,
                " low_budget_min=",lowcheck.minresource,
                " high_budget_min=",highcheck===nothing ? "N/A" :
                                      highcheck.minresource,
                " selected=",selected,
                " (no household consumption floor)")
    end
    return selected
end

function stage3_flexible_prices!(w::TransitionWorkspace{T},
                                 ctx::Stage3Context{T},
                                 x::AbstractVector{T}) where {T}
    inp=ctx.economy
    p=inp.params
    n=length(inp.B)
    length(x)==2*n-1 ||
        throw(DimensionMismatch("flexible-rate transition needs 2T-1 unknowns"))
    w.K[1]=inp.initial.K
    @inbounds for t in 2:n
        w.K[t]=x[t-1]
    end
    @inbounds for t in 1:n
        w.L[t]=x[n-1+t]
        w.K[t]>0 && w.L[t]>0 ||
            throw(DomainError((w.K[t],w.L[t]),"positive capital/labor needed"))
        K=w.K[t]
        w.Y[t]=inp.z[t]*K^p.alpha*w.L[t]^(one(T)-p.alpha)
        mc=(p.epsilon_f-one(T))/p.epsilon_f
        w.rental[t]=mc*p.alpha*w.Y[t]/K
        w.wage[t]=mc*(one(T)-p.alpha)*w.Y[t]/w.L[t]
        w.profits[t]=w.Y[t]/p.epsilon_f
        w.I[t]=t<n ? w.K[t+1]-(one(T)-p.delta)*K : p.delta*K
        if t<n
            gap=w.I[t]/K-p.delta
            w.q[t]=one(T)+inp.phi_ac*gap
            w.adjustment[t]=inp.phi_ac/T(2)*K*gap^2
            w.profits[t]+=(w.q[t]-one(T))*w.I[t]-w.adjustment[t]
        else
            w.q[t]=one(T)
            w.adjustment[t]=zero(T)
        end
        w.q[t]>0 || throw(DomainError(w.q[t],"negative capital shadow price"))
    end
    @inbounds for t in 1:n
        prevq=t==1 ? inp.q_initial : w.q[t-1]
        w.rate[t]=(w.rental[t]+(one(T)-p.delta)*w.q[t])/prevq-one(T)
        one(T)+w.rate[t]>0 ||
            throw(DomainError(w.rate[t],"nonpositive gross return"))
        if t==1 && ctx.foreign && inp.B_initial!=zero(T)
            # Solve the initial trade price with the ACTUAL foreign bond
            # revaluation. Do not first evaluate an unrelated domestic-
            # debt price equation: it may have no feasible positive root.
            w.terms[t]=stage3_positive_trade_root(w,ctx)
            startassets=inp.q_initial*w.K[t]+
                        inp.B_initial*ctx.terms_pre/w.terms[t]
        else
            startassets=t==1 ?
                inp.q_initial*w.K[t]+inp.B_initial :
                inp.B[t-1]+w.q[t-1]*w.K[t]
        end
        endassets=t<n ? inp.B[t]+w.q[t]*w.K[t+1] :
                        inp.B[t]+inp.terminal.K
        w.home_budget[t]=p.chi*
            (startassets*(one(T)+w.rate[t])-endassets+
             w.wage[t]*w.L[t]+w.profits[t])
        available=if t<n
            w.Y[t]-w.I[t]-w.adjustment[t]-w.home_budget[t]
        else
            terminalgap=w.I[t-1]/w.K[t]-p.delta
            cost=inp.phi_ac/T(2)*w.K[t]*terminalgap^2
            w.Y[t]-w.I[t-1]-cost-w.home_budget[t]
        end
        if t==1 && ctx.foreign && inp.B_initial!=zero(T)
            price_error=available-inp.reference.export_shift*
                        w.terms[t]^(-p.theta_star)
            abs(price_error)<T(1e-8) ||
                throw(DomainError(price_error,"initial FX trade root mismatch"))
        else
            available>0 ||
                throw(DomainError(available,"no positive trade-price root"))
            w.terms[t]=(inp.reference.export_shift/available)^
                       (inv(p.theta_star))
        end
        w.home_supply[t]=w.home_budget[t]
        isfinite(w.terms[t]) && w.terms[t]>0 ||
            throw(DomainError(w.terms[t],"invalid terms of trade"))
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
    # Truncated infinite-horizon GE:
    # Dates 1:T-1 are economic transition periods. Date T is the
    # exogenous stationary continuation POLICY boundary. Counting:
    # (T-1) asset + (T-1) labor + (T-1) goods + two terminal
    # anchors (L_T=L_ss and ToT_T=ToT_ss) = 3T-1 equations.
    # DO NOT impose period-T goods/labor clearing against a terminal
    # stationary policy evaluated at still-transitory period-T prices.
    @inbounds for t in 1:n-1
        c=w.aggregate_home[t]
        c>0 || throw(DomainError(c,"nonpositive home consumption"))
        # Paper eq. (5): W/P_H = (epsilon_w/(epsilon_w-1)) *
        # chi_dis * L^varphi * C_H / chi (gamma=theta=1).
        # Taking L = ratio is valid ONLY for varphi=1. The author's
        # MATLAB transition omits the 1/varphi power although its
        # stationary union condition correctly includes L^varphi.
        ratio=unionfactor*w.wage[t]/(inp.reference.chi_dis*c)
        ratio>zero(T) || throw(DomainError(ratio,"nonpositive union labor FOC"))
        w.labor_implied[t]=ratio^(inv(p.varphi))
        res[n-1+t]=(w.labor_implied[t]-w.L[t])/scaleL
        if fixed
            inflation=t==1 ? w.terms[1]/ctx.terms_pre :
                              w.terms[t]/w.terms[t-1]
            nominal_cost=ctx.theta_nr/T(2)*(inflation-one(T))^2
            goods=w.Y[t]-w.I[t]-w.adjustment[t]-nominal_cost-
                  inp.reference.export_shift*w.terms[t]^(-p.theta_star)
            res[2n-1+t]=(goods-c)/max(one(T),inp.reference.Y)
        end
    end
    w.labor_implied[n]=w.L[n]
    res[2n-1]=(w.L[n]-inp.terminal.L)/scaleL
    if fixed
        res[3n-1]=(w.terms[n]-inp.terminal.terms)/
                   max(one(T),abs(inp.terminal.terms))
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

# Numerical Stage 3 solver: finite-difference Jacobian, Broyden updates,
# and adaptive globalization for the same 399/599 GE equations.
# Persist accepted but NOT YET CONVERGED GE paths, separately from
# verified economic outputs. The next run reevaluates all residuals.
function stage3_write_checkpoint(path::AbstractString,
                                 x::AbstractVector{T},
                                 residual::AbstractVector{T},
                                 iteration::Int) where {T<:AbstractFloat}
    length(x)==length(residual) && length(x)>0 ||
        throw(DimensionMismatch("checkpoint needs paired GE vectors"))
    iteration>0 || throw(ArgumentError("checkpoint iteration must be positive"))
    all(isfinite,x) && all(isfinite,residual) ||
        throw(ArgumentError("checkpoint contains nonfinite values"))
    temp=string(path,".tmp")
    open(temp,"w") do io
        println(io,"# STAGE3_UNCONVERGED_CHECKPOINT v1")
        println(io,"# n=",length(x))
        println(io,"# iter=",iteration)
        println(io,"# max_scaled_residual=",maximum(abs,residual))
        for value in x
            println(io,repr(value))
        end
    end
    mv(temp,path;force=true)
    return path
end

function stage3_read_checkpoint(path::AbstractString,n::Int)
    n>0 || throw(ArgumentError("checkpoint dimension must be positive"))
    lines=readlines(path)
    length(lines)==n+4 ||
        throw(DimensionMismatch("checkpoint dimension differs from current GE"))
    lines[1]=="# STAGE3_UNCONVERGED_CHECKPOINT v1" ||
        throw(ArgumentError("invalid Stage 3 checkpoint signature"))
    lines[2]=="# n=$(n)" ||
        throw(DimensionMismatch("checkpoint source horizon differs"))
    startswith(lines[3],"# iter=") &&
    startswith(lines[4],"# max_scaled_residual=") ||
        throw(ArgumentError("checkpoint metadata invalid"))
    x=parse.(Float64,lines[5:end])
    all(isfinite,x) || throw(ArgumentError("nonfinite checkpoint vector"))
    return x
end

"""
A Levenberg–Marquardt direction for the SAME GE equations. This is only
a globalization fallback after a freshly recomputed Jacobian fails to
produce a decreasing Newton step. A nonzero residual is never accepted
as an equilibrium.
"""
function stage3_lm_direction(J::Matrix{T}, f::Vector{T}, factor::T) where {T}
    n=length(f)
    size(J)==(n,n) || throw(DimensionMismatch("square GE Jacobian"))
    factor>zero(T) || throw(ArgumentError("LM damping must be positive"))
    gram=transpose(J)*J
    scale=max(one(T),sum(abs2,J)/T(n))
    @inbounds for j in 1:n
        gram[j,j]+=factor*scale
    end
    rhs=-(transpose(J)*f)
    direction=gram\rhs
    all(isfinite,direction) ||
        throw(ErrorException("nonfinite regularized GE direction"))
    return direction
end

function stage3_nonlinear_solve(ctx::Stage3Context{T},x0::Vector{T};
                                fixed::Bool=false,explicit_flexible::Bool=false,
                                tolerance::T=T(1e-6),
                                maxiter::Int=70,refresh::Int=15,
                                max_evaluations::Int=5000,
                                on_accepted_step::Union{Nothing,Function}=nothing,
                                verbose::Bool=true) where {T}
    n=length(x0)
    fixed && explicit_flexible && throw(ArgumentError("conflicting GE regimes"))
    expected=(fixed || explicit_flexible) ? 3*length(ctx.economy.B)-1 :
                     2*length(ctx.economy.B)-1
    n==expected || throw(DimensionMismatch("Stage3 unknown dimensions"))
    refresh>0 || throw(ArgumentError("Jacobian refresh period must be positive"))
    max_evaluations>0 ||
        throw(ArgumentError("max_evaluations must be positive"))
    work=TransitionWorkspace(ctx.economy)
    # Original MATLAB MIT_transition.m uses MaxFunEvals=5000 for
    # the unexpected credit contraction. Count every trial and every
    # Jacobian probe; never silently continue tens of thousands of
    # expensive household GE evaluations.
    evaluation! = function (res,x)
        work.function_evals<max_evaluations ||
            throw(ErrorException("Stage3 MATLAB-aligned residual evaluation " *
                "budget exhausted at $(work.function_evals)/$(max_evaluations), " *
                "before 399/599-equation convergence"))
        if explicit_flexible
            stage3_flexible_full_residual!(res,x,ctx,work)
        elseif fixed
            stage3_fixed_residual!(res,x,ctx,work)
        else
            stage3_residual!(res,x,ctx,work)
        end
        return nothing
    end
    x=copy(x0)
    f=zeros(T,n)
    trialres=similar(f)
    trial=similar(x)
    jacprobe=similar(x)
    J=Matrix{T}(undef,n,n)
    evaluation!(f,x)
    small_progress=0

    # The original MATLAB fsolve uses a numerical Jacobian and
    # globalization. A Broyden-only Jacobian can become uninformative
    # after repeated small damped steps in the piecewise-linear wealth
    # mass-splitting problem. Refresh *on observed stagnation*.
    function rebuild_jacobian!()
        @inbounds for j in 1:n
            copyto!(jacprobe,x)
            dx=T(1e-5)*max(one(T),abs(x[j]))
            jacprobe[j]+=dx
            valid=true
            try
                evaluation!(trialres,jacprobe)
            catch err
                if err isa DomainError || err isa ArgumentError
                    valid=false
                else
                    rethrow()
                end
            end
            if !valid
                jacprobe[j]=x[j]-dx
                try
                    evaluation!(trialres,jacprobe)
                catch err
                    if err isa DomainError || err isa ArgumentError
                        throw(ErrorException(
                            "Jacobian column $j lacks a feasible finite perturbation; " *
                            "current GE candidate lies on a feasibility boundary"))
                    end
                    rethrow()
                end
                dx=-dx
            end
            @inbounds for i in 1:n
                J[i,j]=(trialres[i]-f[i])/dx
            end
        end
        # EGM and distributions are stateful; restore the *center*
        # evaluation rather than leaving the last Jacobian probe in work.
        evaluation!(f,x)
        return nothing
    end

    function try_step!(direction::Vector{T})
        all(isfinite,direction) || return (false,zero(T),0)
        base=norm(f)
        α=one(T)
        rejected=0
        for _ in 1:22
            @inbounds for j in 1:n
                trial[j]=x[j]+α*direction[j]
            end
            valid=true
            try
                evaluation!(trialres,trial)
            catch err
                if err isa DomainError || err isa ArgumentError
                    valid=false
                else
                    rethrow()
                end
            end
            if valid && all(isfinite,trialres) &&
               norm(trialres)<(one(T)-T(1e-4)*α)*base
                return (true,α,rejected)
            end
            rejected+=1
            α/=T(2)
        end
        return (false,zero(T),rejected)
    end

    for iter in 1:maxiter
        maxres=maximum(abs,f)
        maxat=argmax(abs.(f))
        tag=maxat<length(ctx.economy.B) ? "asset" :
            maxat<2*length(ctx.economy.B) ? "labor" : "goods"
        oldnorm=norm(f)
        verbose && println("STAGE3 ",explicit_flexible ? "FLEX_FULL" :
                           fixed ? "FIXED" : "FLEX",
                           " iteration=",iter-1,
                           " scaled_max=",maxres," max_at=",maxat,
                           " market=",tag," norm=",oldnorm,
                           " evaluations=",work.function_evals)
        if maxres<=tolerance
            evaluation!(f,x)
            maximum(abs,f)<=tolerance ||
                error("Stage3 root residual changed on repeat evaluation")
            return TransitionSolveResult(ctx.economy,work,copy(x),copy(f),
                        iter-1,work.function_evals,maximum(abs,f),true)
        end
        did_refresh=false
        # A tiny accepted alpha is NOT evidence that an expensive 399-column
        # Jacobian is stale. Use measured residual progress, not alpha alone,
        # to decide when another full matrix must be computed.
        if iter==1 || (iter-1)%refresh==0 || small_progress>=2
            verbose && println("STAGE3 JACOBIAN refresh at iteration=",iter-1,
                               " reason=",iter==1 ? "initial" :
                               small_progress>=2 ? "stagnation" : "scheduled")
            rebuild_jacobian!()
            did_refresh=true
            small_progress=0
        end
        direction=-(J\f)
        accepted,α,rejects=try_step!(direction)

        if !accepted && !did_refresh
            # Refresh a STALE Broyden Jacobian at most once in this
            # iteration. If it was just computed, try LM directly.
            verbose && println("STAGE3 Newton line search failed; " *
                               "recomputing full numerical Jacobian")
            rebuild_jacobian!()
            did_refresh=true
            direction=-(J\f)
            accepted,α,rejects=try_step!(direction)
        end

        if !accepted
            for λ in (T(1e-6),T(1e-4),T(1e-2),T(1e-1),one(T),T(10),T(100))
                direction=stage3_lm_direction(J,f,λ)
                accepted,α,rejects=try_step!(direction)
                if accepted
                    verbose && println("STAGE3 LM globalization accepted λ=",λ,
                                       " damping=",α)
                    break
                end
            end
        end

        accepted || throw(ErrorException(
            "Stage3 cannot reduce full GE residual even with fresh Jacobian " *
            "and source-compatible regularized steps. " *
            "max_index=$maxat max_level=$maxres norm=$oldnorm " *
            "evaluations=$(work.function_evals). " *
            "Do not treat this as a converged equilibrium."))
        dx=α.*direction
        change=trialres.-f
        correction=change.-J*dx
        denominator=dot(dx,dx)
        denominator>zero(T) || error("vanishing Stage3 GE step")
        @inbounds for j in 1:n
            weight=dx[j]/denominator
            for i in 1:n
                J[i,j]+=correction[i]*weight
            end
        end
        newnorm=norm(trialres)
        progress=one(T)-newnorm/oldnorm
        small_progress=progress<T(0.02) ? small_progress+1 : 0
        verbose && println("STAGE3 STEP alpha=",α," rejected=",rejects,
                           " relative_norm_improvement=",progress,
                           " next_max=",maximum(abs,trialres))
        copyto!(x,trial)
        copyto!(f,trialres)
        # The checkpoint is explicitly NOT a converged equilibrium. Only
        # fully accepted iterations are persisted, after f and x agree.
        if on_accepted_step!==nothing
            on_accepted_step(x,f,iter)
        end
    end
    evaluation!(f,x)
    throw(ErrorException(
        "Stage3 GE not converged after $maxiter iterations: " *
        "max=$(maximum(abs,f)), evaluations=$(work.function_evals)"))
end

function solve_stage3_flexible(ctx::Stage3Context{T};
                               x0::AbstractVector{T}=initial_transition_guess(ctx.economy),
                               kwargs...) where {T}
    return stage3_nonlinear_solve(ctx,Vector{T}(x0);fixed=false,kwargs...)
end
