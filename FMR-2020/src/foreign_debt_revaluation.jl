
"""
Source-level implementation of MIT_transition.m and backsolve_egm.m:
inherited capital, pre-shock nominal capital price, gross household capital
and foreign liabilities, and unexpected exchange-rate balance-sheet revaluation.
No household net-wealth sign is silently discarded or truncated.
"""
struct Stage3Context{T<:AbstractFloat}
    economy::TransitionInputs{T}
    mu_pre::Matrix{T}
    terms_pre::T
    leverage_hat::T
    theta_nr::T
    foreign::Bool
    label::String
end

function stage3_full_grids(gr::AssetGrids{T}) where {T}
    # Authors' Agrid_fine includes pre-revaluation negative wealth points.
    # The stationary grids are unchanged; this view is only for MIT shocks.
    return AssetGrids(gr.coarse_full,gr.fine_full,gr.coarse,
                      gr.fine_full,gr.borrowing_limit,gr.maximum_asset)
end

function stage3_expand_distribution(mu::AbstractMatrix{T},
                                    original::Vector{T},
                                    all_assets::Vector{T}) where {T}
    size(mu,1)==length(original) ||
        throw(DimensionMismatch("pre-shock distribution/asset knots"))
    ans=zeros(T,length(all_assets),size(mu,2))
    @inbounds for i in eachindex(original)
        k=searchsortedfirst(all_assets,original[i])
        k<=length(all_assets) && all_assets[k]==original[i] ||
            throw(ArgumentError("coarse and full fine asset knots differ"))
        for j in axes(mu,2)
            ans[k,j]=mu[i,j]
        end
    end
    abs(sum(ans)-sum(mu))<T(1e-11) ||
        error("expanding the initial distribution loses mass")
    return ans
end

function stage3_initial_state(st::SteadyState{T},K::T,L::T,B::T,
                              qprev::T,termspre::T,
                              mu::Matrix{T}) where {T}
    return SteadyState(st.name,B,st.beta,st.r,K,L,st.Y,st.wage,termspre,
                       st.investment,qprev,st.rental,st.profits,
                       st.consumption_home,st.consumption_foreign,
                       st.aggregate_assets,st.asset_market_error,
                       st.goods_market_error,st.union_error,
                       st.household,mu,st.distribution_iterations)
end

function stage3_expansion_context(steady;theta_nr::Float64=100.0)
    base=transition_inputs(steady)
    gr=stage3_full_grids(base.grids)
    mu=stage3_expand_distribution(steady.initial.distribution,
                                  base.grids.fine,gr.fine)
    initial=stage3_initial_state(steady.initial,steady.initial.K,
                                 steady.initial.L,steady.initial.B,
                                 one(steady.initial.K),steady.initial.terms,mu)
    inp=TransitionInputs(base.params,base.reference,base.chain,gr,
                         initial,base.terminal,base.B,base.z,base.phi_ac,
                         one(steady.initial.K),steady.initial.B)
    return Stage3Context(inp,copy(mu),steady.initial.terms,0.0,
                         theta_nr,true,"foreign-debt expansion")
end

"""
    credit_reversal_path(previous, Bcal; shock_start=41, periods=200,
                         decay_periods=periods)

Return a credit reversal that reaches Bcal at date decay_periods.
When extending the numerical transition horizon, hold decay_periods
fixed so the economic shock at ALL original dates is unchanged.
Increasing periods while also increasing decay_periods would define
a different economic experiment, not a valid horizon-convergence test.
"""
function credit_reversal_path(previous::AbstractVector{T}, Bcal::T;
                              shock_start::Int=41, periods::Int=200,
                              decay_periods::Int=periods) where {T}
    2<=shock_start<=length(previous) ||
        throw(ArgumentError("shock must be after first observed quarter"))
    periods>1 && decay_periods>1 ||
        throw(ArgumentError("reversal and decay horizons must exceed one"))
    # MATLAB shock_start_t=41; shock_b_agg(shock_start_t-1).
    bprev=previous[shock_start-1]
    return [Bcal+(bprev-Bcal)*
            max(zero(T),one(T)-T(j-1)/T(decay_periods-1))^4
            for j in 1:periods]
end

"""
    stage3_extend_horizon(ctx, periods)

Lengthen an EXISTING Stage 3 economic experiment, preserving the
entire original foreign-credit path, shocked wealth distribution,
inherited balance sheet and policy regime. The additional dates
continue at the limiting exogenous B and productivity z. This is a
horizon-truncation experiment, NOT a recalibration of the shock.
"""
function stage3_extend_horizon(ctx::Stage3Context{T}, periods::Int) where {T}
    old=ctx.economy
    previous=length(old.B)
    periods>=previous ||
        throw(ArgumentError("new Stage3 horizon cannot be shorter"))
    if periods==previous
        return ctx
    end
    # The original 200-date credit expansion is calibrated against
    # stationary B = B_200, not the AR(1) infinite limit. Preserving the
    # economic experiment requires maintaining that terminal supply after
    # date 200; extending the AR(1) tail while keeping terminal SS fixed
    # would violate terminal asset market conditions.
    abs(old.B[end]-old.terminal.B)<=T(1e-8)*
        max(one(T),abs(old.terminal.B)) ||
        throw(ArgumentError("cannot append stationary B when source shock has not reached its calibrated terminal"))
    B=vcat(old.B,fill(old.terminal.B,periods-previous))
    z=vcat(old.z,fill(old.z[end],periods-previous))
    extended=TransitionInputs(old.params,old.reference,old.chain,
                              old.grids,old.initial,old.terminal,B,z,
                              old.phi_ac,old.q_initial,old.B_initial)
    return Stage3Context(extended,copy(ctx.mu_pre),ctx.terms_pre,
                         ctx.leverage_hat,ctx.theta_nr,ctx.foreign,
                         ctx.label*" | horizon="*string(periods))
end

"""
    stage3_prolong_aggregate_guess(old_path, old_ctx, extended_ctx; rho=0.9)

Prolong a 3T-1 K/L/terms numerical starting vector to a longer horizon
without changing any originally solved period. New capital dates approach
the unchanged terminal stationary capital geometrically; newly appended
L and terms dates start at their stationary values. This changes ONLY
the numerical initial guess, never the exogenous experiment or equations.
"""
function stage3_prolong_aggregate_guess(
    old_path::AbstractVector{T}, old_ctx::Stage3Context{T},
    extended_ctx::Stage3Context{T};
    rho::T=T(0.9)) where {T<:AbstractFloat}

    old_n=length(old_ctx.economy.B)
    new_n=length(extended_ctx.economy.B)
    length(old_path)==3old_n-1 ||
        throw(DimensionMismatch("source aggregate path must be 3T-1"))
    new_n>=old_n ||
        throw(ArgumentError("extended horizon must not be shorter"))
    zero(T)<rho<one(T) ||
        throw(ArgumentError("prolongation rho must lie in (0,1)"))
    extended_ctx.economy.B[1:old_n]==old_ctx.economy.B ||
        throw(ArgumentError("prolongation must preserve the credit experiment"))
    K0=old_ctx.economy.initial.K
    K0==extended_ctx.economy.initial.K ||
        throw(ArgumentError("initial capital cannot change on extension"))
    Kold=vcat(K0,old_path[1:old_n-1])
    Lold=old_path[old_n:2old_n-1]
    xold=old_path[2old_n:3old_n-1]
    ntail=new_n-old_n
    Ktail=[extended_ctx.economy.terminal.K+
           (Kold[end]-extended_ctx.economy.terminal.K)*rho^s
           for s in 1:ntail]
    Ltail=fill(extended_ctx.economy.terminal.L,ntail)
    xtail=fill(extended_ctx.economy.terminal.terms,ntail)
    ans=vcat(Kold[2:end],Ktail,Lold,Ltail,xold,xtail)
    length(ans)==3new_n-1 ||
        throw(DimensionMismatch("prolonged 3T-1 unknown vector"))
    all(x->isfinite(x)&&x>zero(T),ans) ||
        throw(DomainError(ans,"nonpositive prolonged aggregate candidate"))
    return ans
end

function stage3_reversal_context(steady,previous::TransitionSolveResult{T};
                                shock_start::Int=41,
                                leverage_hat::T=T(1/16),
                                theta_nr::T=T(100)) where {T}
    previous.converged || throw(ArgumentError("pre-shock GE must converge"))
    w=previous.workspace
    2<=shock_start<=length(previous.input.B) ||
        throw(ArgumentError("invalid reversal date"))
    allgr=stage3_full_grids(steady.grids)
    grid_previous=previous.input.grids.fine
    pre=grid_previous==allgr.fine ?
        copy(view(w.distribution,:,:,shock_start)) :
        stage3_expand_distribution(view(w.distribution,:,:,shock_start),
                                   grid_previous,allgr.fine)
    bpre=previous.input.B[shock_start-1]
    qpre=w.q[shock_start-1]
    xpre=w.terms[shock_start-1]
    Bnew=credit_reversal_path(previous.input.B,steady.calibrated.B;
                              shock_start=shock_start,
                              periods=steady.parameters.transition_horizon)
    initial=stage3_initial_state(steady.initial,w.K[shock_start],
                                 w.L[shock_start],bpre,qpre,xpre,pre)
    inp=TransitionInputs(steady.parameters,steady.reference,steady.chain,
                         allgr,initial,steady.calibrated,Bnew,
                         ones(T,length(Bnew)),T(17),qpre,bpre)
    return Stage3Context(inp,pre,xpre,leverage_hat,theta_nr,true,
                         "unexpected reversal")
end

"""
Backsolve_egm.m capital/bond mapping, including source special case BG_shock=0.
Revalued net wealth is capital + (bond / (x_new/x_pre)).
The Jacobian of the wealth mapping is NOT a probability Jacobian:
discrete mass is moved via neighboring-node weights, preserving total mass.

The original gross-capital allocation does NOT generally satisfy
sum(mu_i*k_i) == q_previous*K_initial when leverage_hat > 0. The explicit
opt-in keyword strict_portfolio=true rescales ALL gross capital positions
by the common factor (q_previous*K_initial)/sum(mu_i*k_i) BEFORE computing
foreign bonds b_i=a_i-k_i. This preserves their relative capital exposure,
ensures aggregate gross capital equals q_previous*K_initial, and makes
total bonds equal the inherited B when the pre-shock asset market clears.
The default remains false to reproduce and audit the earlier saved roots.
Those roots are NOT equilibria of the strict-portfolio system.
"""
function revalue_foreign_debt!(destination::AbstractMatrix{T},
                               original::AbstractMatrix{T},
                               grid::AbstractVector{T},
                               Kinit::T, qprev::T, Bprev::T,
                               xprev::T, xnew::T, leverage_hat::T;
                               foreign::Bool=true,
                               strict_portfolio::Bool=false) where {T}
    size(destination)==size(original) ||
        throw(DimensionMismatch("revaluation distribution"))
    size(original,1)==length(grid) ||
        throw(DimensionMismatch("wealth grid"))
    xprev>0 && xnew>0 || throw(DomainError((xprev,xnew),"positive exchange rates"))
    0<=leverage_hat<1 ||
        throw(ArgumentError("leverage_hat must be in [0,1)"))
    sum(original)>=0 || throw(ArgumentError("negative probability"))
    original===destination &&
        throw(ArgumentError("source/destination alias in revaluation"))
    rich=zero(T)
    @inbounds for j in axes(original,2),i in eachindex(grid)
        if grid[i]>0
            rich+=original[i,j]*grid[i]
        end
    end
    lev=Bprev==zero(T) ? one(T) :
        (one(T)-leverage_hat)*Kinit*qprev/rich
    isfinite(lev) && lev>0 ||
        throw(DomainError(lev,"invalid gross-asset leverage"))
    kfixed=leverage_hat*Kinit
    target_capital=qprev*Kinit
    target_capital>zero(T) ||
        throw(DomainError(target_capital,"positive inherited capital value required"))
    # The MATLAB-style raw allocation can under/over-allocate aggregate
    # capital due to its fixed intercept and max(.,0) truncation.
    # A single common scalar keeps relative portfolio exposures intact.
    raw_capital=zero(T)
    pre_wealth=zero(T)
    @inbounds for j in axes(original,2),i in eachindex(grid)
        mass=original[i,j]
        mass==zero(T) && continue
        raw_capital+=mass*max(lev*grid[i]+kfixed,zero(T))
        pre_wealth+=mass*grid[i]
    end
    raw_capital>zero(T) && isfinite(raw_capital) ||
        throw(DomainError(raw_capital,"gross capital allocation cannot be zero"))
    capital_scale=strict_portfolio ? target_capital/raw_capital : one(T)
    if strict_portfolio
        # This checks the INHERITED t-1 asset market, not the future
        # market. Only tiny numerical equilibrium error is admissible.
        mismatch=pre_wealth-(target_capital+Bprev)
        abs(mismatch)<=T(1e-5)*max(one(T),abs(pre_wealth),
                                  abs(target_capital+Bprev)) ||
            throw(DomainError(mismatch,
                "pre-shock wealth does not equal inherited qK+B"))
    end
    ratio=xnew/xprev
    fill!(destination,zero(T))
    total_capital=zero(T)
    total_bonds=zero(T)
    total_revalued=zero(T)
    @inbounds for j in axes(original,2),i in eachindex(grid)
        mass=original[i,j]
        mass==0 && continue
        mass>0 || throw(DomainError(mass,"negative pre-shock probability"))
        a=grid[i]
        cap=capital_scale*max(lev*a+kfixed,zero(T))
        bonds=a-cap
        new_a=foreign ? cap+bonds/ratio : a
        if new_a<grid[1]-T(1e-9) || new_a>grid[end]+T(1e-9)
            throw(DomainError(new_a,
              "positive-mass foreign-debt revaluation leaves full asset grid"))
        end
        if new_a<=grid[1]
            idx=1
            weight=one(T)
        elseif new_a>=grid[end]
            idx=length(grid)-1
            weight=zero(T)
        else
            idx=min(length(grid)-1,searchsortedlast(grid,new_a))
            weight=(grid[idx+1]-new_a)/(grid[idx+1]-grid[idx])
        end
        destination[idx,j]+=weight*mass
        destination[idx+1,j]+=(one(T)-weight)*mass
        total_capital+=mass*cap
        total_bonds+=mass*bonds
        total_revalued+=mass*new_a
    end
    before=sum(original)
    after=sum(destination)
    abs(before-after)<=T(1e-10) ||
        throw(ErrorException("foreign revaluation loses household probability"))
    actual_wealth=zero(T)
    @inbounds for j in axes(destination,2),i in eachindex(grid)
        actual_wealth+=destination[i,j]*grid[i]
    end
    abs(actual_wealth-total_revalued)<=T(1e-8)*
         max(one(T),abs(total_revalued)) ||
        throw(ErrorException("wealth revaluation inconsistent with mass interpolation"))
    if strict_portfolio
        abs(total_capital-target_capital)<=T(1e-9)*
             max(one(T),abs(target_capital)) ||
            throw(ErrorException("strict gross-capital market does not clear"))
        abs(total_bonds-Bprev)<=T(1e-5)*
             max(one(T),abs(Bprev)) ||
            throw(ErrorException("strict inherited foreign bond market does not clear"))
    end
    return (mass_before=before,mass_after=after,
            capital=total_capital,bonds=total_bonds,
            capital_target=target_capital,
            capital_market_gap=total_capital-target_capital,
            inherited_bond_target=Bprev,
            inherited_bond_market_gap=total_bonds-Bprev,
            legacy_gross_capital_gap=raw_capital-target_capital,
            capital_normalization=capital_scale,
            strict_portfolio=strict_portfolio,
            net_wealth_before=total_capital+total_bonds,
            net_wealth_after=actual_wealth,
            predicted_net_wealth=total_capital+total_bonds/ratio,
            exchange_ratio=ratio,leverage=lev,fixed_capital=kfixed)
end
