# Reference translation of author solve_R.m and backsolve_egm.m.
# This is deliberately NOT substituted into stage3_household_backward!
# until the MATLAB and Julia one-step numerical arrays are compared.
# The MATLAB numerical consumption floor is recorded as a source feature,
# not assumed to satisfy the household's economic budget constraint.

function source_matlab_grid_layout(gr::AssetGrids{T}) where {T<:AbstractFloat}
    b=gr.borrowing_limit
    function borrow_index(grid)
        i=argmin(abs.(grid .- b))
        abs(grid[i]-b)<=T(1e-10)*max(one(T),abs(b)) ||
            throw(ArgumentError("source asset grid lacks borrowing-limit node"))
        return i
    end
    ic=borrow_index(gr.coarse_full)
    jf=borrow_index(gr.fine_full)
    return (coarse_total=length(gr.coarse_full),
            fine_total=length(gr.fine_full),
            coarse_borrow_index=ic,
            fine_borrow_index=jf,
            coarse_below=ic-1,
            fine_below=jf-1,
            coarse_active=length(gr.coarse_full)-ic+1,
            fine_active=length(gr.fine_full)-jf+1)
end

# Exactly one backward-EGM date with the ORIGINAL full coarse/fine asset
# grids and Index_b_min, using MATLAB's source-consumption-floor convention.
# Inputs are the full-grid date-(t+1) consumption, date-t firm prices,
# date-t and t+1 real interest rates, and date-t terms of trade.
# This one-period reference map is separate from the budget-strict map.
function source_matlab_egm_step(gr::AssetGrids{T},
                                chain::ProductivityChain{T},
                                p::ModelParams{T},
                                future::AbstractMatrix{T},
                                beta::T, rate_now::T, rate_next::T,
                                terms::T, wage::T, labor::T, profits::T;
                                consumption_floor::T=T(1e-5)) where {T<:AbstractFloat}
    coarse=gr.coarse_full
    fine=gr.fine_full
    layout=source_matlab_grid_layout(gr)
    nc=length(coarse)
    nf=length(fine)
    ns=length(chain.levels)
    size(future)==(nc,ns) ||
        throw(DimensionMismatch("MATLAB backsolve_egm requires full 193-row terminal cpol"))
    all(isfinite,future) && minimum(future)>zero(T) ||
        throw(ArgumentError("source future consumption must be positive"))
    beta>zero(T) && one(T)+rate_now>zero(T) &&
        one(T)+rate_next>zero(T) && consumption_floor>zero(T) ||
        throw(ArgumentError("source EGM parameters invalid"))
    m=consumption_cost(terms,p)
    m>zero(T) || throw(DomainError(m,"nonpositive source expenditure factor"))
    now=Matrix{T}(undef,nc,ns)
    cf=Matrix{T}(undef,nf,ns)
    af=Matrix{T}(undef,nf,ns)
    raw_budget=Matrix{T}(undef,nf,ns)
    cutoff=Vector{T}(undef,ns)
    a_today=Vector{T}(undef,nc)
    c_s=similar(a_today)
    slopes=similar(a_today)
    constrained=similar(a_today)
    constrained_slopes=similar(a_today)
    gross_now=one(T)+rate_now
    gross_next=one(T)+rate_next
    for j in 1:ns
        income=wage*labor*chain.levels[j]+profits
        for i in 1:nc
            mu=zero(T)
            for jp in 1:ns
                mu+=chain.transition[j,jp]*future[i,jp]^(-p.gamma)
            end
            c_s[i]=(beta*gross_next*mu)^(-inv(p.gamma))
            a_today[i]=(m*c_s[i]+coarse[i]-income)/gross_now
            constrained[i]=max((gross_now*coarse[i]+income-
                               gr.borrowing_limit)/m,consumption_floor)
        end
        all(diff(a_today).>zero(T)) ||
            throw(ArgumentError("nonmonotone full-source endogenous asset nodes"))
        pchip_slopes!(slopes,a_today,c_s)
        pchip_slopes!(constrained_slopes,coarse,constrained)
        threshold=a_today[layout.coarse_borrow_index]
        cutoff[j]=threshold
        for i in 1:nc
            a=coarse[i]
            now[i,j]=a<=threshold ? constrained[i] :
                max(pchip_evaluate(a_today,c_s,slopes,a),consumption_floor)
        end
        for i in 1:nf
            a=fine[i]
            binding=a<=threshold
            c=binding ?
                pchip_evaluate(coarse,constrained,constrained_slopes,a) :
                pchip_evaluate(a_today,c_s,slopes,a)
            cf[i,j]=c
            af[i,j]=binding ? gr.borrowing_limit :
                gross_now*a+income-m*c
            raw_budget[i,j]=gross_now*a+income-m*c-af[i,j]
        end
    end
    return (coarse_consumption=now,fine_consumption=cf,
            fine_savings=af,budget_gap=raw_budget,
            endogenous_threshold=cutoff,layout=layout)
end

# Recover exactly the unscaled MATLAB market equations from the
# existing Julia normalized 399-/599-vector. This does NOT imply the
# economic residual maps are equivalent: full-grid EGM and branch
# selection still differ and require source-based comparison.
function stage3_unscale_residual!(dest::AbstractVector{T},
                                  normalized::AbstractVector{T},
                                  ctx::Stage3Context{T};
                                  fixed::Bool=false) where {T<:AbstractFloat}
    n=length(ctx.economy.B)
    dim=fixed ? 3n-1 : 2n-1
    length(dest)==dim && length(normalized)==dim ||
        throw(DimensionMismatch("399/599 market residual dimension mismatch"))
    scaleK=max(one(T),ctx.economy.initial.K)
    scaleL=max(one(T),ctx.economy.initial.L)
    @inbounds for t in 1:n-1
        dest[t]=normalized[t]*scaleK
    end
    @inbounds for t in n:(2n-1)
        dest[t]=normalized[t]*scaleL
    end
    if fixed
        scaleY=max(one(T),ctx.economy.reference.Y)
        @inbounds for t in 2n:3n-1
            dest[t]=normalized[t]*scaleY
        end
    end
    return dest
end
