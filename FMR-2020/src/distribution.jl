"""
Matrix-free household distribution transition. Current and next arrays are
asset x productivity. Each incumbent mass is split between adjacent asset
knots, then over idiosyncratic productivity transitions.
"""
function distribution_step!(next::AbstractMatrix{T}, current::AbstractMatrix{T},
                            savings::AbstractMatrix{T}, grid::AbstractVector{T},
                            P::AbstractMatrix{T}) where {T<:AbstractFloat}
    n,ns = size(current)
    size(next) == (n,ns) || throw(DimensionMismatch("next distribution"))
    size(savings) == (n,ns) || throw(DimensionMismatch("savings policy"))
    size(P) == (ns,ns) || throw(DimensionMismatch("productivity transition"))
    length(grid) == n || throw(DimensionMismatch("asset grid"))
    fill!(next,zero(T))
    @inbounds for j in 1:ns, i in 1:n
        mass = current[i,j]
        mass == 0 && continue
        mass >= 0 || throw(DomainError(mass,"negative mass"))
        save = savings[i,j]
        if save < grid[1]-T(1e-8) || save > grid[end]+T(1e-8)
            throw(DomainError(save,"savings outside active distribution grid"))
        end
        if save <= grid[1]
            low = 1
            weight = one(T)
        elseif save >= grid[end]
            low = n-1
            weight = zero(T)
        else
            low = min(n-1,searchsortedlast(grid,save))
            weight = (grid[low+1]-save)/(grid[low+1]-grid[low])
        end
        if !(zero(T) <= weight <= one(T))
            throw(DomainError(weight,"invalid asset probability"))
        end
        for jp in 1:ns
            increment = mass*P[j,jp]
            next[low,jp] += weight*increment
            next[low+1,jp] += (one(T)-weight)*increment
        end
    end
    return next
end

"""
Precomputed interpolation brackets and weights for fixed household saving
policies. Reused during invariant-distribution power iterations.
"""
struct AssetProjection{T<:AbstractFloat}
    lower::Matrix{Int}
    weight::Matrix{T}
end

function asset_projection(savings::AbstractMatrix{T},
                          grid::AbstractVector{T}) where {T<:AbstractFloat}
    n,ns=size(savings)
    length(grid)==n || throw(DimensionMismatch("grid/policy"))
    lower=Matrix{Int}(undef,n,ns)
    weight=Matrix{T}(undef,n,ns)
    @inbounds for j in 1:ns, i in 1:n
        v=savings[i,j]
        if v < grid[1]-T(1e-8) || v > grid[end]+T(1e-8)
            throw(DomainError(v,"saving leaves the active grid"))
        end
        if v<=grid[1]
            lower[i,j]=1
            weight[i,j]=one(T)
        elseif v>=grid[end]
            lower[i,j]=n-1
            weight[i,j]=zero(T)
        else
            k=min(n-1,searchsortedlast(grid,v))
            lower[i,j]=k
            weight[i,j]=(grid[k+1]-v)/(grid[k+1]-grid[k])
        end
        0 <= weight[i,j] <= 1 || throw(DomainError(weight[i,j],"bad interpolation weight"))
    end
    return AssetProjection(lower,weight)
end

function distribution_step!(next::AbstractMatrix{T},current::AbstractMatrix{T},
                            projected::AssetProjection{T},
                            P::AbstractMatrix{T}) where {T<:AbstractFloat}
    n,ns=size(current)
    size(next)==(n,ns) || throw(DimensionMismatch("distribution"))
    size(projected.lower)==(n,ns) || throw(DimensionMismatch("projection"))
    size(P)==(ns,ns) || throw(DimensionMismatch("Markov matrix"))
    fill!(next,zero(T))
    @inbounds for j in 1:ns, i in 1:n
        mass=current[i,j]
        mass==0 && continue
        mass >= 0 || throw(DomainError(mass,"negative distribution mass"))
        k=projected.lower[i,j]
        w=projected.weight[i,j]
        for jp in 1:ns
            moved=mass*P[j,jp]
            next[k,jp] += w*moved
            next[k+1,jp] += (1-w)*moved
        end
    end
    return next
end

function stationary_distribution(p::ModelParams{T},gr::AssetGrids{T},
                                 chain::ProductivityChain{T},
                                 hh::HouseholdSolution{T}) where {T}
    n,ns = size(hh.savings_fine)
    current = Matrix{T}(undef,n,ns)
    @inbounds for j in 1:ns, i in 1:n
        current[i,j] = chain.invariant[j]/n
    end
    next = similar(current)
    projected = asset_projection(hh.savings_fine,gr.fine)
    last_error = T(Inf)
    for iteration in 1:p.dist_maxiter
        distribution_step!(next,current,projected,chain.transition)
        last_error = maximum(abs.(next .- current))
        current,next = next,current
        if last_error < p.dist_tol
            mass = sum(current)
            abs(mass-one(T)) < T(1e-10) ||
                error("stationary distribution loses probability mass")
            return (distribution=current, iterations=iteration, residual=last_error)
        end
    end
    throw(ErrorException("stationary distribution did not converge; residual $(last_error)"))
end

function distribution_invariance(mu::Matrix{T},hh::HouseholdSolution{T},
                                 gr::AssetGrids{T},chain::ProductivityChain{T}) where {T}
    next = similar(mu)
    distribution_step!(next,mu,hh.savings_fine,gr.fine,chain.transition)
    return (mass=sum(mu), min_mass=minimum(mu),
            invariance_error=maximum(abs.(next .- mu)),
            next_mass=sum(next))
end
