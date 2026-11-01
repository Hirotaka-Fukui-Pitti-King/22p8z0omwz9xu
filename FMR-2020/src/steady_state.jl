struct SteadyState{T<:AbstractFloat}
    name::String
    B::T
    beta::T
    r::T
    K::T
    L::T
    Y::T
    wage::T
    terms::T
    investment::T
    capital_price::T
    rental::T
    profits::T
    consumption_home::T
    consumption_foreign::T
    aggregate_assets::T
    asset_market_error::T
    goods_market_error::T
    union_error::T
    household::HouseholdSolution{T}
    distribution::Matrix{T}
    distribution_iterations::Int
end

"""
Bisection with an explicit sign-changing bracket. A small function update
is not sufficient: the residual is independently checked at the root.
"""
function bracketed_root(f, lo::T, hi::T;
                        tol::T=T(1e-8), maxiter::Int=80) where {T<:AbstractFloat}
    lo < hi || throw(ArgumentError("root bracket must be ordered"))
    flo = f(lo)
    fhi = f(hi)
    isfinite(flo) && isfinite(fhi) ||
        throw(ArgumentError("nonfinite root residual at bracket endpoint"))
    abs(flo) <= tol && return (root=lo, residual=flo, iterations=0)
    abs(fhi) <= tol && return (root=hi, residual=fhi, iterations=0)
    signbit(flo) != signbit(fhi) ||
        throw(ArgumentError("root not bracketed: f(lo)=$(flo), f(hi)=$(fhi)"))
    for iteration in 1:maxiter
        mid = lo+(hi-lo)/2
        fm = f(mid)
        isfinite(fm) || throw(DomainError(fm,"nonfinite residual"))
        if abs(fm) <= tol
            return (root=mid,residual=fm,iterations=iteration)
        end
        if signbit(fm) == signbit(flo)
            lo = mid
            flo = fm
        else
            hi = mid
            fhi = fm
        end
    end
    throw(ErrorException("root solve did not meet tolerance after $(maxiter) iterations"))
end

"""
Given r and B, compute (K,L,Y,w,profits,x) using firm FOCs, the
labor-union condition and the aggregate budget / foreign-demand identity.
This is the steady-state aggregate closure prior to household aggregation.
"""
function aggregate_prices(p::ModelParams{T},ref::ReferenceCalibration{T},
                          r::T,B::T) where {T}
    r+p.delta > 0 || throw(DomainError(r,"r+delta must be positive"))
    markup_inv = (p.epsilon_f-1)/p.epsilon_f
    k_l = ((r+p.delta)/(markup_inv*p.alpha))^(inv(p.alpha-1))
    wage = markup_inv*(1-p.alpha)*k_l^p.alpha
    wm = (p.epsilon_w-1)/p.epsilon_w

    function union_gap(L::T)
        K = k_l*L
        Y = k_l^p.alpha*L
        profits = Y/p.epsilon_f
        ch_budget = p.chi*(r*(K+B)+wage*L+profits)
        ch_union = p.chi*wm*wage/(ref.chi_dis*L^p.varphi)
        return ch_union-ch_budget
    end
    ls = bracketed_root(union_gap,T(1e-5),T(100);tol=p.root_tol,maxiter=p.root_maxiter)
    L = ls.root
    K = k_l*L
    Y = k_l^p.alpha*L
    profits = Y/p.epsilon_f
    ch_budget = p.chi*(r*(K+B)+wage*L+profits)
    ch_budget > 0 || throw(DomainError(ch_budget,"nonpositive household budget consumption"))
    exports = Y-p.delta*K-ch_budget
    exports > 0 || throw(DomainError(exports,"negative steady-state exports"))
    terms = (exports/ref.export_shift)^(-inv(p.theta_star))
    terms > 0 || throw(DomainError(terms,"nonpositive terms of trade"))
    return (K=K,L=L,Y=Y,wage=wage,terms=terms,profits=profits,
            ch_budget=ch_budget,rental=r+p.delta,union_error=ls.residual)
end

function stationary_candidate(p::ModelParams{T},ref::ReferenceCalibration{T},
                              chain::ProductivityChain{T},gr::AssetGrids{T},
                              beta::T,r::T,B::T;
                              name::String="candidate") where {T}
    aggregate = aggregate_prices(p,ref,r,B)
    prices = HouseholdPrices(r,aggregate.wage,aggregate.terms,
                             aggregate.L,aggregate.profits)
    hh = solve_household(p,gr,chain,prices,beta)
    mu_result = stationary_distribution(p,gr,chain,hh)
    mu = mu_result.distribution
    assets = sum(mu .* hh.savings_fine)
    ch = sum(mu .* hh.consumption_fine)
    cf = ((one(T)-p.chi)/p.chi)*aggregate.terms^p.theta*ch
    goods = aggregate.Y - p.delta*aggregate.K -
            ref.export_shift*aggregate.terms^(-p.theta_star) - ch
    asset = assets - (aggregate.K+B)
    return SteadyState(name,B,beta,r,aggregate.K,aggregate.L,aggregate.Y,
                       aggregate.wage,aggregate.terms,p.delta*aggregate.K,
                       one(T),aggregate.rental,aggregate.profits,ch,cf,assets,
                       asset,goods,aggregate.union_error,
                       hh,mu,mu_result.iterations)
end

function steady_state_diagnostics(p::ModelParams{T},chain::ProductivityChain{T},
                                  gr::AssetGrids{T},state::SteadyState{T}) where {T}
    pr = HouseholdPrices(state.r,state.wage,state.terms,state.L,state.profits)
    hd = household_diagnostics(p,gr,chain,pr,state.beta,state.household)
    dd = distribution_invariance(state.distribution,state.household,gr,chain)
    macro_scale = max(one(T),abs(state.K+state.B),abs(state.Y))
    factor_price_error = max(abs(state.rental - (p.epsilon_f-1)/p.epsilon_f*
                                         p.alpha*state.Y/state.K),
                             abs(state.wage - (p.epsilon_f-1)/p.epsilon_f*
                                         (1-p.alpha)*state.Y/state.L))
    return (asset_error=state.asset_market_error,
            scaled_asset_error=abs(state.asset_market_error)/macro_scale,
            goods_error=state.goods_market_error,
            scaled_goods_error=abs(state.goods_market_error)/macro_scale,
            union_error=state.union_error,
            factor_price_error=factor_price_error,
            household=hd,distribution=dd)
end
