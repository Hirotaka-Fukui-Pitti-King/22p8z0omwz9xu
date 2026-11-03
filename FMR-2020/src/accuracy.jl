"""
An independent Euler/KKT diagnostic evaluated at every *fine-grid* household
state, using the stationary fine-grid consumption policy as the continuation
policy. Unlike the existing coarse-grid maximum, this also records the
stationary household probability mass attached to each residual.
"""
struct EulerRecord{T<:AbstractFloat}
    asset_index::Int
    productivity_index::Int
    asset::T
    productivity::T
    stationary_mass::T
    next_asset::T
    consumption::T
    relative_euler_gap::T
    euler_kkt_violation::T
    binding::Bool
end

function euler_accuracy_map(p::ModelParams{T}, gr::AssetGrids{T},
                            chain::ProductivityChain{T},
                            st::SteadyState{T}) where {T<:AbstractFloat}
    cf = st.household.consumption_fine
    sf = st.household.savings_fine
    mu = st.distribution
    n, ns = size(cf)
    size(sf) == (n,ns) && size(mu) == (n,ns) ||
        throw(DimensionMismatch("fine policy/distribution arrays"))
    length(gr.fine) == n && length(chain.levels) == ns ||
        throw(DimensionMismatch("grid/productivity array lengths"))
    slopes = similar(cf)
    @inbounds for j in 1:ns
        pchip_slopes!(view(slopes,:,j),gr.fine,view(cf,:,j))
    end
    records = Vector{EulerRecord{T}}(undef,n*ns)
    biggest_nonbinding = zero(T)
    biggest_kkt = zero(T)
    weighted_sum = zero(T)
    weighted_nonbinding = zero(T)
    positive_mass = zero(T)
    total_mass = zero(T)
    top_asset_tail_mass = zero(T)
    maximum_mass_at_single_state = zero(T)
    threshold = gr.borrowing_limit+T(1e-8)*max(one(T),abs(gr.borrowing_limit))
    idx = 0
    @inbounds for j in 1:ns, i in 1:n
        idx += 1
        current_c = cf[i,j]
        chosen_a = sf[i,j]
        current_c > 0 || throw(DomainError(current_c,"nonpositive fine consumption"))
        expected = zero(T)
        for jp in 1:ns
            future_c = pchip_evaluate(gr.fine,view(cf,:,jp),
                                      view(slopes,:,jp),chosen_a)
            future_c > 0 || throw(DomainError(future_c,"nonpositive future consumption"))
            expected += chain.transition[j,jp]*future_c^(-p.gamma)
        end
        lhs = current_c^(-p.gamma)
        rhs = st.beta*(1+st.r)*expected
        gap = (lhs-rhs)/lhs
        binding = chosen_a <= threshold
        violation = binding ? max(zero(T),-gap) : abs(gap)
        mass = mu[i,j]
        mass >= 0 || throw(DomainError(mass,"negative stationary probability"))
        records[idx] = EulerRecord(i,j,gr.fine[i],chain.levels[j],
                                   mass,chosen_a,current_c,gap,violation,binding)
        weighted_sum += mass*violation
        total_mass += mass
        if mass > 0
            positive_mass += mass
            maximum_mass_at_single_state = max(maximum_mass_at_single_state,mass)
        end
        if binding
            biggest_kkt = max(biggest_kkt,violation)
        else
            biggest_nonbinding = max(biggest_nonbinding,violation)
            weighted_nonbinding += mass*violation
        end
        if gr.fine[i] >= T(0.9)*gr.fine[end]
            top_asset_tail_mass += mass
        end
    end
    abs(total_mass-one(T)) <= T(1e-8) ||
        throw(ArgumentError("Euler diagnostics require unit-mass distribution"))
    ranked = sort(records;by=x->x.euler_kkt_violation,rev=true)
    mass_order = sort(records;by=x->x.euler_kkt_violation)
    accumulated = zero(T)
    p99 = T(NaN)
    for record in mass_order
        accumulated += record.stationary_mass
        if accumulated >= T(0.99)*total_mass
            p99 = record.euler_kkt_violation
            break
        end
    end
    summary = (max_unconstrained=biggest_nonbinding,
               max_kkt_violation=biggest_kkt,
               probability_weighted_mean=weighted_sum/total_mass,
               probability_weighted_unconstrained=weighted_nonbinding/total_mass,
               probability_weighted_p99=p99,
               mass_in_top_ten_percent_asset_grid=top_asset_tail_mass,
               max_state_probability=maximum_mass_at_single_state,
               total_mass=total_mass,
               n_active_states=count(x->x.stationary_mass > 0,records),
               n_grid_states=length(records))
    return (summary=summary,ranked=ranked,records=records)
end


"""
    household_best_response_saving(grid, future_consumption, probabilities,
                                   gamma, beta, gross_next, resources,
                                   unit_cost, borrowing_limit; atol=1e-11)

Solve a single household's one-period constrained saving problem conditional
on independently supplied next-date consumption policies. The current
budget is unit_cost*c + saving = resources. For the implemented paper
calibration (gamma=theta=1), the interior Euler condition uses
c^(-gamma) = beta*gross_next*E[c_next^(-gamma)].

Unlike inverse-EGM interpolation, this brackets the Euler root directly in
the economically feasible interval. The borrowing limit binds when the
Euler surplus at that limit is nonnegative. It never clips a choice to an
artificial upper asset-grid bound: if the optimum lies beyond the supported
future-policy grid, it throws instead. Intended as an independent, targeted
accuracy check, not a new aggregate GE solver.
"""
function household_best_response_saving(
    grid::AbstractVector{T}, future_consumption::AbstractMatrix{T},
    probabilities::AbstractVector{T}, gamma::T, beta::T, gross_next::T,
    resources::T, unit_cost::T, borrowing_limit::T;
    atol::T=T(1e-11)) where {T<:AbstractFloat}

    na=length(grid)
    ns=length(probabilities)
    na>=2 && size(future_consumption)==(na,ns) ||
        throw(DimensionMismatch("future household policy/grid mismatch"))
    all(diff(grid).>zero(T)) ||
        throw(ArgumentError("future-policy grid must be increasing"))
    abs(sum(probabilities)-one(T))<=T(1e-10) &&
        all(x->isfinite(x)&&x>=zero(T),probabilities) ||
        throw(ArgumentError("productivity transition probabilities must sum to one"))
    gamma>zero(T) && beta>zero(T) && gross_next>zero(T) &&
        unit_cost>zero(T) && atol>zero(T) ||
        throw(ArgumentError("invalid household Euler parameters"))
    all(x->isfinite(x)&&x>zero(T),future_consumption) ||
        throw(DomainError(future_consumption,
                          "future consumption must be finite and positive"))
    grid[1]<=borrowing_limit<=grid[end] ||
        throw(DomainError(borrowing_limit,"borrowing limit outside grid"))

    # Numerical interior bracket, not an economic consumption floor.
    # A true policy solution requires positive consumption and a' >= a_min.
    maximum_saving=min(grid[end],
                       resources-unit_cost*T(1e-12))
    maximum_saving>borrowing_limit ||
        throw(DomainError(resources,"no positive feasible consumption"))

    slopes=Matrix{T}(undef,na,ns)
    for jp in 1:ns
        pchip_slopes!(view(slopes,:,jp),grid,
                      view(future_consumption,:,jp))
    end
    function surplus(saving::T)
        c=(resources-saving)/unit_cost
        c>zero(T) ||
            throw(DomainError(c,"nonpositive household consumption"))
        emu=zero(T)
        for jp in 1:ns
            pnext=pchip_evaluate(grid,view(future_consumption,:,jp),
                                 view(slopes,:,jp),saving)
            pnext>zero(T) ||
                throw(DomainError(pnext,"nonpositive interpolated continuation"))
            emu+=probabilities[jp]*pnext^(-gamma)
        end
        return c^(-gamma)-beta*gross_next*emu
    end

    low_gap=surplus(borrowing_limit)
    if low_gap>=zero(T)
        c=(resources-borrowing_limit)/unit_cost
        return (savings=borrowing_limit,consumption=c,
                binding=true,euler_surplus=low_gap,iterations=0)
    end

    upper_gap=surplus(maximum_saving)
    upper_gap>zero(T) ||
        throw(DomainError((maximum_saving,upper_gap),
                          "unconstrained optimal saving exceeds the interpolated future-policy grid"))
    lo=borrowing_limit
    hi=maximum_saving
    iterations=0
    for it in 1:100
        mid=lo+(hi-lo)/T(2)
        val=surplus(mid)
        if val>zero(T)
            hi=mid
        else
            lo=mid
        end
        iterations=it
        if hi-lo<=atol*max(one(T),abs(mid))
            saving=lo+(hi-lo)/T(2)
            return (savings=saving,
                    consumption=(resources-saving)/unit_cost,
                    binding=false,euler_surplus=surplus(saving),
                    iterations=iterations)
        end
    end
    throw(ErrorException("independent household Euler root did not converge"))
end


"""
    euler_refine_one_saving(grid, future, slopes, probs, gamma, beta,
                            gross_next, resources, m, amin, guess;
                            trigger=1e-5, target=1e-9)

A safeguarded one-dimensional Euler/KKT corrector for a *given* next-date
consumption function. Uses PRECOMPUTED PCHIP slopes, so unlike the
single-state audit's standalone solver it avoids allocating slopes at
every asset state. Returns (saving, consumption, refined, relative_gap,
offgrid). It never clamps a genuinely off-grid optimum: offgrid=true
leaves the original candidate unchanged for subsequent feasibility
checking. This routine is a fixed-price household numerical refinement,
not an aggregate equilibrium solver.
"""
function euler_refine_one_saving(
    grid::AbstractVector{T},future::AbstractMatrix{T},
    slopes::AbstractMatrix{T},probs::AbstractVector{T},
    gamma::T,beta::T,gross_next::T,resources::T,m::T,
    amin::T,guess::T;trigger::T=T(1e-5),
    target::T=T(1e-9)) where {T<:AbstractFloat}

    ns=length(probs)
    size(future)==size(slopes)==(length(grid),ns) ||
        throw(DimensionMismatch("future-policy interpolation dimensions"))
    gamma>zero(T) && beta>zero(T) && gross_next>zero(T) &&
        m>zero(T) && trigger>target>zero(T) ||
        throw(ArgumentError("invalid Euler corrector parameters"))
    lo=amin
    hi=min(grid[end],resources-m*T(1e-12))
    if !isfinite(guess) || guess<lo-T(1e-9)
        throw(DomainError(guess,"infeasible pre-correction saving"))
    end
    if hi<=lo
        throw(DomainError(resources,"no positive constrained consumption"))
    end
    if guess>hi || guess<grid[1]
        c=(resources-guess)/m
        return (saving=guess,consumption=c,refined=false,
                relative_gap=T(NaN),offgrid=true)
    end
    euler(s)=begin
        c=(resources-s)/m
        c>zero(T) || throw(DomainError(c,"nonpositive Euler trial consumption"))
        emu=zero(T)
        @inbounds for jp in 1:ns
            fc=pchip_evaluate(grid,view(future,:,jp),
                              view(slopes,:,jp),s)
            fc>zero(T) || throw(DomainError(fc,"nonpositive continuation consumption"))
            emu+=probs[jp]*fc^(-gamma)
        end
        return c^(-gamma)-beta*gross_next*emu
    end
    c0=(resources-guess)/m
    g=euler(guess)
    relative=g/(c0^(-gamma))
    if abs(relative)<=trigger ||
       (guess<=lo+T(1e-9) && relative>=-trigger)
        return (saving=guess,consumption=c0,refined=false,
                relative_gap=relative,offgrid=false)
    end

    glo=euler(lo)
    if glo>=zero(T)
        c=(resources-lo)/m
        return (saving=lo,consumption=c,refined=true,
                relative_gap=glo/c^(-gamma),offgrid=false)
    end
    ghi=euler(hi)
    if ghi<=zero(T)
        # No root on the supported continuation-policy grid.
        # Never manufacture a solution by clipping to grid[end].
        return (saving=guess,consumption=c0,refined=false,
                relative_gap=relative,offgrid=true)
    end
    saving=clamp(guess,lo,hi)
    for _ in 1:32
        c=(resources-saving)/m
        g=euler(saving)
        relative=g/c^(-gamma)
        if abs(relative)<=target
            return (saving=saving,consumption=c,refined=true,
                    relative_gap=relative,offgrid=false)
        end
        if g>zero(T)
            hi=saving
        else
            lo=saving
        end
        step=T(1e-5)*max(one(T),abs(saving))
        a=max(lo,saving-step)
        b=min(hi,saving+step)
        proposal=(lo+hi)/T(2)
        if b>a
            derivative=(euler(b)-euler(a))/(b-a)
            if isfinite(derivative) && derivative>zero(T)
                newton=saving-g/derivative
                if isfinite(newton) && lo<newton<hi
                    proposal=newton
                end
            end
        end
        saving=proposal
    end
    throw(ErrorException("Euler/KKT local Newton-bisection corrector did not converge"))
end
