
"""
T periods correspond to author MATLAB j=1,...,TT.
K[1] is inherited capital; K[2:T] are unknown. L[1:T] unknown.
B[t] is end-of-period foreign-asset supply, not opening-period B.
"""
struct TransitionInputs{T<:AbstractFloat}
    params::ModelParams{T}
    reference::ReferenceCalibration{T}
    chain::ProductivityChain{T}
    grids::AssetGrids{T}
    initial::SteadyState{T}
    terminal::SteadyState{T}
    B::Vector{T}
    z::Vector{T}
    phi_ac::T
    q_initial::T
    B_initial::T
end

function credit_expansion_path(p::ModelParams{T}, Y_reference::T;
                               periods::Int=p.transition_horizon) where {T}
    periods > 1 || throw(ArgumentError("transition requires T>=2"))
    path = Vector{T}(undef,periods)
    debt = zero(T)
    @inbounds for t in 1:periods
        debt -= p.current_account_initial*Y_reference*p.credit_persistence^(t-1)
        path[t] = debt
    end
    return path
end

function transition_inputs(steady;
                           periods::Int=steady.parameters.transition_horizon,
                           phi_ac::Float64=17.0, zero_shock::Bool=false)
    p=steady.parameters
    T=typeof(steady.initial.K)
    B=zero_shock ? fill(steady.initial.B,periods) :
                   credit_expansion_path(p,steady.reference.Y;periods=periods)
    terminal=zero_shock ? steady.initial : steady.final
    if !zero_shock && periods != p.transition_horizon
        throw(ArgumentError("new horizon requires a newly calibrated terminal steady state"))
    end
    return TransitionInputs(p,steady.reference,steady.chain,steady.grids,
                            steady.initial,terminal,B,ones(T,periods),
                            T(phi_ac),one(T),steady.initial.B)
end

mutable struct TransitionWorkspace{T<:AbstractFloat}
    K::Vector{T}
    L::Vector{T}
    Y::Vector{T}
    I::Vector{T}
    adjustment::Vector{T}
    q::Vector{T}
    rental::Vector{T}
    wage::Vector{T}
    rate::Vector{T}
    profits::Vector{T}
    terms::Vector{T}
    home_budget::Vector{T}
    home_supply::Vector{T}
    aggregate_home::Vector{T}
    aggregate_savings::Vector{T}
    labor_implied::Vector{T}
    coarse_consumption::Array{T,3}
    fine_consumption::Array{T,3}
    fine_savings::Array{T,3}
    distribution::Array{T,3}
    endogenous_asset::Vector{T}
    endogenous_consumption::Vector{T}
    interpolation_derivative::Vector{T}
    distribution_error::Vector{T}
    function_evals::Int
end

function TransitionWorkspace(input::TransitionInputs{T}) where {T}
    n=length(input.B)
    nc=length(input.grids.coarse)
    nf=length(input.grids.fine)
    ns=length(input.chain.levels)
    a()=zeros(T,n)
    return TransitionWorkspace(a(),a(),a(),a(),a(),a(),a(),a(),a(),
                               a(),a(),a(),a(),a(),a(),a(),
                               zeros(T,nc,ns,n),zeros(T,nf,ns,n),
                               zeros(T,nf,ns,n),zeros(T,nf,ns,n),
                               zeros(T,nc),zeros(T,nc),zeros(T,nc),
                               a(),0)
end

function initial_transition_guess(input::TransitionInputs{T}) where {T}
    n=length(input.B)
    x=Vector{T}(undef,2n-1)
    @inbounds for t in 2:n
        fraction=one(T)-T(0.9)^(t-1)
        x[t-1]=input.initial.K+fraction*(input.terminal.K-input.initial.K)
    end
    @inbounds for t in 1:n
        fraction=one(T)-T(0.9)^(t-1)
        x[n-1+t]=input.initial.L+fraction*(input.terminal.L-input.initial.L)
    end
    return x
end

"""
Compute flexible-price aggregate paths. The author solves one scalar
terms-of-trade equation per date. With theta=1, m(x)=1/chi, so this
equation has an exact positive scalar root; using it is algebraically
identical to a converged scalar fsolve.
"""
function transition_prices!(w::TransitionWorkspace{T},
                            input::TransitionInputs{T},
                            x::AbstractVector{T}) where {T}
    p=input.params
    n=length(input.B)
    length(x)==2n-1 || throw(DimensionMismatch("expected 2T-1 unknowns"))
    w.K[1]=input.initial.K
    @inbounds for t in 2:n
        w.K[t]=x[t-1]
    end
    @inbounds for t in 1:n
        w.L[t]=x[n-1+t]
        w.K[t]>0 && w.L[t]>0 ||
            throw(DomainError((w.K[t],w.L[t]),"K,L must be positive"))
    end
    @inbounds for t in 1:n
        K=w.K[t]
        L=w.L[t]
        Y=input.z[t]*K^p.alpha*L^(1-p.alpha)
        w.Y[t]=Y
        w.rental[t]=(p.epsilon_f-1)/p.epsilon_f*p.alpha*Y/K
        w.wage[t]=(p.epsilon_f-1)/p.epsilon_f*(1-p.alpha)*Y/L
        w.profits[t]=Y/p.epsilon_f
        w.I[t]=t<n ? w.K[t+1]-(1-p.delta)*K : p.delta*K
        if t<n
            gap=w.I[t]/K-p.delta
            w.adjustment[t]=input.phi_ac/2*K*gap^2
            w.q[t]=1+input.phi_ac*gap
            w.profits[t]+=(w.q[t]-1)*w.I[t]-w.adjustment[t]
        else
            w.adjustment[t]=zero(T)
            w.q[t]=one(T)
        end
        w.q[t]>0 || throw(DomainError(w.q[t],"capital price must be positive"))
    end
    @inbounds for t in 1:n
        previous_q=t==1 ? input.q_initial : w.q[t-1]
        w.rate[t]=(w.rental[t]+(1-p.delta)*w.q[t])/previous_q-1
        1+w.rate[t]>0 || throw(DomainError(w.rate[t],"gross return nonpositive"))
        beginning=t==1 ? input.q_initial*w.K[1]+input.B_initial :
                          input.B[t-1]+w.q[t-1]*w.K[t]
        ending=t==n ? input.B[t]+input.terminal.K :
                       input.B[t]+w.q[t]*w.K[t+1]
        w.home_budget[t]=p.chi*(beginning*(1+w.rate[t])-ending+
                                w.wage[t]*w.L[t]+w.profits[t])
        available=if t<n
            w.Y[t]-w.I[t]-w.adjustment[t]-w.home_budget[t]
        else
            gap=w.I[t-1]/w.K[t]-p.delta
            cost=input.phi_ac/2*w.K[t]*gap^2
            w.Y[t]-w.I[t-1]-cost-w.home_budget[t]
        end
        available>0 || throw(DomainError(available,"no positive trade-price root"))
        w.terms[t]=(input.reference.export_shift/available)^(inv(p.theta_star))
        w.home_supply[t]=w.home_budget[t]
        isfinite(w.terms[t]) && w.terms[t]>0 ||
            throw(DomainError(w.terms[t],"invalid terms of trade"))
    end
    return w
end
