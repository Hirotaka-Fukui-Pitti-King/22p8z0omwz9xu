struct AssetGrids{T<:AbstractFloat}
    coarse_full::Vector{T}
    fine_full::Vector{T}
    coarse::Vector{T}
    fine::Vector{T}
    borrowing_limit::T
    maximum_asset::T
end

function make_asset_grids(p::ModelParams{T}, K_reference::T) where {T}
    K_reference > 0 || throw(ArgumentError("positive reference capital required"))
    lower = -p.borrowing_limit_scale*K_reference
    upper = p.asset_max_scale*K_reference
    function build(n::Int)
        z = range(p.grid_min_z,one(T),length=n)
        full = [lower+(upper-lower)*v^p.grid_power for v in z]
        all(diff(full) .> 0) || error("asset grid must be strictly increasing")
        any(isapprox(a,lower;atol=1e-10,rtol=0) for a in full) ||
            throw(ArgumentError("grid does not contain the borrowing limit"))
        return full
    end
    coarse_full = build(p.n_coarse)
    fine_full = build(p.n_fine)
    # Only feasible post-decision assets enter stationary household policies.
    # Negative pre-revaluation states are preserved in *_full for future MIT stages.
    coarse = [a for a in coarse_full if a >= lower-1e-10]
    fine = [a for a in fine_full if a >= lower-1e-10]
    coarse[1] = lower
    fine[1] = lower
    return AssetGrids(coarse_full,fine_full,coarse,fine,lower,upper)
end
