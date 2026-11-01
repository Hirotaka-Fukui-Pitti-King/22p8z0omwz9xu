struct HouseholdPrices{T<:AbstractFloat}
    r::T
    wage::T
    terms::T
    labor::T
    profits::T
end

struct HouseholdSolution{T<:AbstractFloat}
    consumption_coarse::Matrix{T}
    savings_coarse::Matrix{T}
    consumption_fine::Matrix{T}
    savings_fine::Matrix{T}
    iterations::Int
    update_error::T
end

@inline consumption_cost(x::T, p::ModelParams{T}) where {T} =
    one(T)+(one(T)-p.chi)/p.chi * x^(p.theta-one(T))

@inline labor_income(pr::HouseholdPrices{T}, productivity::T) where {T} =
    pr.wage*pr.labor*productivity

"""
Solve the stationary household problem, with active assets a >= borrowing limit.
Full grids (including below-limit MIT-shock states) remain in AssetGrids for Stage 2.
Shape-preserving Hermite interpolation is used on endogenous asset knots.
"""
function solve_household(p::ModelParams{T}, gr::AssetGrids{T},
                         chain::ProductivityChain{T}, pr::HouseholdPrices{T},
                         beta::T) where {T<:AbstractFloat}
    0 < beta < 1 || throw(ArgumentError("beta must lie in (0,1)"))
    pr.r > -1 || throw(ArgumentError("gross real return must be positive"))
    pr.terms > 0 && pr.labor > 0 || throw(ArgumentError("terms and labor must be positive"))
    n = length(gr.coarse)
    ns = length(chain.levels)
    length(gr.fine) >= n || throw(ArgumentError("fine grid smaller than coarse grid"))
    m = consumption_cost(pr.terms,p)
    m > 0 || throw(ArgumentError("nonpositive consumption price"))
    gross = 1+pr.r
    c = Matrix{T}(undef,n,ns)
    fresh = similar(c)
    @inbounds for j in 1:ns, i in 1:n
        resources = gross*gr.coarse[i] + labor_income(pr,chain.levels[j]) + pr.profits
        cap = (resources-gr.borrowing_limit)/m
        cap > 0 || throw(DomainError(cap,"stationary state without feasible consumption"))
        c[i,j] = min(max((T(0.03)*gr.coarse[i]+labor_income(pr,chain.levels[j])+pr.profits)/m,
                          T(1e-12)),cap)
    end
    endo_asset = Vector{T}(undef,n)
    endo_consumption = similar(endo_asset)
    derivatives = similar(endo_asset)
    last_error = T(Inf)
    completed = false
    n_iter = 0

    for iter in 1:p.egm_maxiter
        for j in 1:ns
            income = labor_income(pr,chain.levels[j]) + pr.profits
            @inbounds for i in 1:n
                emu = zero(T)
                for jp in 1:ns
                    emu += chain.transition[j,jp] * c[i,jp]^(-p.gamma)
                end
                endo_consumption[i] = (beta*gross*emu)^(-inv(p.gamma))
                endo_asset[i] = (m*endo_consumption[i]+gr.coarse[i]-income)/gross
            end
            all(diff(endo_asset) .> 0) ||
                throw(ArgumentError("EGM endogenous grid not increasing"))
            pchip_slopes!(derivatives,endo_asset,endo_consumption)
            @inbounds for i in 1:n
                a = gr.coarse[i]
                new_c = if a <= endo_asset[1]
                    (gross*a+income-gr.borrowing_limit)/m
                else
                    pchip_evaluate(endo_asset,endo_consumption,derivatives,a)
                end
                isfinite(new_c) && new_c > 0 ||
                    throw(DomainError(new_c,"invalid EGM consumption"))
                fresh[i,j] = new_c
            end
        end
        last_error = maximum(abs.(fresh .- c))
        n_iter = iter
        if last_error <= p.egm_tol
            copyto!(c,fresh)
            completed = true
            break
        end
        weight = last_error < T(1e-4) ? T(0.5) : T(0.9)
        @inbounds for j in 1:ns, i in 1:n
            c[i,j] = weight*fresh[i,j] + (one(T)-weight)*c[i,j]
        end
    end
    completed || throw(ErrorException("EGM failed to converge after $(n_iter) iterations; residual $(last_error)"))
    savings = similar(c)
    @inbounds for j in 1:ns, i in 1:n
        a = gr.coarse[i]
        sav = gross*a + labor_income(pr,chain.levels[j]) + pr.profits - m*c[i,j]
        if sav < gr.borrowing_limit-T(1e-8)*max(one(T),abs(gr.borrowing_limit))
            throw(DomainError(sav,"computed savings violate the borrowing limit"))
        end
        savings[i,j] = max(sav,gr.borrowing_limit)
        if savings[i,j] == gr.borrowing_limit
            c[i,j] = (gross*a+labor_income(pr,chain.levels[j])+pr.profits-gr.borrowing_limit)/m
        end
    end
    nf = length(gr.fine)
    cf = Matrix{T}(undef,nf,ns)
    gf = similar(cf)

    # Author's stationary EGM: evaluate fine-grid consumption DIRECTLY
    # from the endogenous current-asset and consumption knots, as in
    # solve_R.m and solve_EGM_EL_open_tg.m (not a second interpolation
    # of a coarse-grid saving policy). Keep Julia-native array kernels.
    #
    # At a <= endogenous asset knot for post-decision a_min, the borrowing
    # constraint binds and consumption comes directly from the budget.
    # Nonbinding consumption is PCHIP-interpolated on endogenous knots.
    # All saving decisions follow from the same household budget identity.
    @inbounds for j in 1:ns
        income = labor_income(pr,chain.levels[j]) + pr.profits
        for i in 1:n
            expected_mu = zero(T)
            for jp in 1:ns
                expected_mu += chain.transition[j,jp]*c[i,jp]^(-p.gamma)
            end
            endo_consumption[i] = (beta*gross*expected_mu)^(-inv(p.gamma))
            endo_asset[i] = (m*endo_consumption[i]+gr.coarse[i]-income)/gross
        end
        all(diff(endo_asset) .> 0) ||
            throw(ArgumentError("EGM endogenous asset grid is not increasing"))
        pchip_slopes!(derivatives,endo_asset,endo_consumption)
        for i in 1:nf
            a = gr.fine[i]
            if a <= endo_asset[1]
                sav = gr.borrowing_limit
                consum = (gross*a+income-sav)/m
            else
                consum = pchip_evaluate(endo_asset,endo_consumption,
                                       derivatives,a)
                sav = gross*a+income-m*consum
            end
            isfinite(consum) && consum > zero(T) ||
                throw(DomainError(consum,
                    "EGM direct fine-grid consumption must be positive"))
            isfinite(sav) && sav >= gr.borrowing_limit ||
                throw(DomainError(sav,
                    "EGM direct fine-grid savings below borrowing limit"))
            sav <= gr.fine[end] ||
                throw(DomainError(sav,
                    "fine-grid policy exceeds distribution grid upper bound"))
            cf[i,j] = consum
            gf[i,j] = sav
        end
    end
    return HouseholdSolution(c,savings,cf,gf,n_iter,last_error)
end

"""
Check Euler/KKT, monotone savings, and budget identity on the coarse grid.
The Euler expectation is evaluated at each optimal a' by independent Hermite
interpolation of the resulting consumption policy.
"""
function household_diagnostics(p::ModelParams{T}, gr::AssetGrids{T},
                               chain::ProductivityChain{T},pr::HouseholdPrices{T},
                               beta::T,hh::HouseholdSolution{T}) where {T}
    n,ns = size(hh.consumption_coarse)
    derivatives = similar(gr.coarse)
    all_slopes = Matrix{T}(undef,n,ns)
    for j in 1:ns
        pchip_slopes!(derivatives,gr.coarse,view(hh.consumption_coarse,:,j))
        all_slopes[:,j] .= derivatives
    end
    max_euler_unconstrained = zero(T)
    max_kkt_violation = zero(T)
    max_budget = zero(T)
    max_fine_budget = zero(T)
    minimum_fine_saving = T(Inf)
    minimum_savings_step = T(Inf)
    @inbounds for j in 1:ns, i in 1:n
        a = gr.coarse[i]
        aprime = hh.savings_coarse[i,j]
        c = hh.consumption_coarse[i,j]
        lhs = c^(-p.gamma)
        expected = zero(T)
        for jp in 1:ns
            cnext = pchip_evaluate(gr.coarse,view(hh.consumption_coarse,:,jp),
                                   view(all_slopes,:,jp),aprime)
            cnext > 0 || throw(DomainError(cnext,"Euler extrapolation not positive"))
            expected += chain.transition[j,jp]*cnext^(-p.gamma)
        end
        rhs = beta*(1+pr.r)*expected
        gap = (lhs-rhs)/lhs
        if aprime <= gr.borrowing_limit+T(1e-8)*max(one(T),abs(gr.borrowing_limit))
            max_kkt_violation = max(max_kkt_violation,max(-gap,zero(T)))
        else
            max_euler_unconstrained = max(max_euler_unconstrained,abs(gap))
        end
        gap_budget = (1+pr.r)*a+labor_income(pr,chain.levels[j])+pr.profits -
                     consumption_cost(pr.terms,p)*c-aprime
        max_budget = max(max_budget,abs(gap_budget))
        if i > 1
            minimum_savings_step = min(minimum_savings_step,
                                       hh.savings_coarse[i,j]-hh.savings_coarse[i-1,j])
        end
    end
    @inbounds for j in 1:ns, i in eachindex(gr.fine)
        a = gr.fine[i]
        sav = hh.savings_fine[i,j]
        consum = hh.consumption_fine[i,j]
        minimum_fine_saving = min(minimum_fine_saving,sav)
        balance = (1+pr.r)*a+labor_income(pr,chain.levels[j])+pr.profits -
                  consumption_cost(pr.terms,p)*consum-sav
        max_fine_budget = max(max_fine_budget,abs(balance))
    end
    return (max_euler_unconstrained=max_euler_unconstrained,
            max_kkt_violation=max_kkt_violation,
            max_budget_error=max_budget,
            max_fine_budget_error=max_fine_budget,
            minimum_fine_saving=minimum_fine_saving,
            minimum_savings_step=minimum_savings_step,
            policy_iterations=hh.iterations,
            policy_update_error=hh.update_error)
end
