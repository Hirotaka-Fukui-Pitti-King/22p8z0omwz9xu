
"""
Forward transport of the full joint asset-productivity distribution.
Policy functions are recomputed by backward EGM for EVERY outer
aggregate-path candidate. The terminal household policy is fixed only
as a boundary condition, never along the interior transition.
"""
function transition_distribution_forward!(w::TransitionWorkspace{T},
                                          input::TransitionInputs{T}) where {T}
    n=length(input.B)
    copyto!(view(w.distribution,:,:,1),input.initial.distribution)
    nstate=length(input.grids.fine)
    ns=length(input.chain.levels)
    @inbounds for t in 1:n
        mu=view(w.distribution,:,:,t)
        ch=zero(T)
        aend=zero(T)
        for j in 1:ns, i in 1:nstate
            mass=mu[i,j]
            ch+=mass*w.fine_consumption[i,j,t]
            aend+=mass*w.fine_savings[i,j,t]
        end
        w.aggregate_home[t]=ch
        w.aggregate_savings[t]=aend
        w.distribution_error[t]=abs(sum(mu)-one(T))
        if w.distribution_error[t]>T(1e-8) || minimum(mu)<0
            throw(DomainError(w.distribution_error[t],
                              "invalid transition probability distribution at t=$t"))
        end
        if t<n
            distribution_step!(view(w.distribution,:,:,t+1),
                               mu,view(w.fine_savings,:,:,t),
                               input.grids.fine,input.chain.transition)
        end
    end
    return w
end
