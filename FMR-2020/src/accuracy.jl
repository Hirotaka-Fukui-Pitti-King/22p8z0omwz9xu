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
