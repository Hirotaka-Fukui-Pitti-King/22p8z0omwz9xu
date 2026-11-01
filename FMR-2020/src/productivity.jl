struct ProductivityChain{T<:AbstractFloat}
    levels::Vector{T}
    transition::Matrix{T}
    invariant::Vector{T}
end

"""
Rouwenhorst (1995) approximation of an AR(1) in log productivity.
The input sigma_stationary is the unconditional standard deviation in logs.
No external package is needed. Rows of the transition matrix sum to one.
"""
function rouwenhorst(n::Int, rho::T, sigma_stationary::T) where {T<:AbstractFloat}
    n >= 2 || throw(ArgumentError("Rouwenhorst requires n >= 2"))
    abs(rho) < 1 || throw(ArgumentError("abs(rho) must be < 1"))
    sigma_stationary > 0 || throw(ArgumentError("sigma must be positive"))
    up = (one(T)+rho)/2
    down = up
    P = T[up 1-up; 1-down down]
    for m in 3:n
        prev = P
        P = zeros(T,m,m)
        @inbounds for j in 1:(m-1), i in 1:(m-1)
            v = prev[i,j]
            P[i,j] += up*v
            P[i,j+1] += (1-up)*v
            P[i+1,j] += (1-down)*v
            P[i+1,j+1] += down*v
        end
        @views P[2:(m-1),:] ./= 2
    end
    @inbounds for i in 1:n
        abs(sum(@view P[i,:])-1) < 1e-12 ||
            error("Rouwenhorst rows do not sum to one")
    end

    invariant = fill(inv(T(n)),n)
    next = similar(invariant)
    converged = false
    for _ in 1:200000
        fill!(next,zero(T))
        @inbounds for j in 1:n, i in 1:n
            next[j] += invariant[i]*P[i,j]
        end
        residual = maximum(abs.(next .- invariant))
        copyto!(invariant,next)
        if residual < 1e-14
            converged = true
            break
        end
    end
    converged || error("productivity invariant distribution did not converge")
    invariant ./= sum(invariant)
    extent = sqrt(T(n-1))*sigma_stationary
    logs = collect(range(-extent,extent,length=n))
    levels = exp.(logs)
    levels ./= sum(invariant .* levels)
    return ProductivityChain(levels,P,invariant)
end

rouwenhorst(p::ModelParams{T}) where {T} =
    rouwenhorst(p.n_s,p.rho_s,p.sigma_log_stationary)
