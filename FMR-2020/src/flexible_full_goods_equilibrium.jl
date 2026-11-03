# The original flexible sudden-stop solve nests a scalar goods-market
# fsolve inside 399 asset/labor equations. An actual author-saved K/L
# vector has no positive first-date scalar goods root: fsolve returns a
# least-squares minimizer, and the outer 399 residual can nevertheless be
# approximately zero. The explicit formulation instead includes all T
# terms-of-trade variables and all T ACTUAL household goods-market equations.
#
# Unknowns = K[2:T], L[1:T], terms[1:T] (3T-1 = 599 for T=200).
# Equations = (T-1) asset + T labor + T home-goods clearing.
#
# This is an ECONOMICALLY CLOSED, distinct model from the original nested
# 399-equation source-replication algorithm, not a claim that the author's
# stored 399-vector is an economic equilibrium.

"""
    stage3_flexible_full_prices!(w, ctx, x)

Generate flexible-price firm returns, wages, profits and the household
terms-of-trade path for an explicit 3T-1 candidate. This reuses the existing
fixed-regime price routine with the Rotemberg adjustment parameter set to
ZERO. Its NKPC then reduces IDENTICALLY to the flexible markup:
mc = (epsilon_f-1)/epsilon_f at every date. Investment costs, the external
debt valuation and the terminal capital condition are unchanged.
"""
function stage3_flexible_full_prices!(w::TransitionWorkspace{T},
                                      ctx::Stage3Context{T},
                                      x::AbstractVector{T}) where {T<:AbstractFloat}
    n=length(ctx.economy.B)
    length(x)==3n-1 ||
        throw(DimensionMismatch("explicit flexible GE requires 3T-1 unknowns"))
    flexible_context=Stage3Context(ctx.economy,ctx.mu_pre,ctx.terms_pre,
                                   ctx.leverage_hat,zero(T),ctx.foreign,
                                   ctx.label)
    stage3_fixed_prices!(w,flexible_context,x)
    return w
end

"""
    stage3_flexible_full_market_residual!(res, w, ctx)

Use the same 399 household asset/labor equations as existing Stage 3, and
APPEND 200 explicit, realized household goods clearing equations. The
last-date boundary is identical to the established fixed-rate transition,
except there is no nominal adjustment cost under flexible pricing.
"""
function stage3_flexible_full_market_residual!(res::AbstractVector{T},
                                               w::TransitionWorkspace{T},
                                               ctx::Stage3Context{T}) where {T<:AbstractFloat}
    inp=ctx.economy
    p=inp.params
    n=length(inp.B)
    length(res)==3n-1 ||
        throw(DimensionMismatch("explicit flexible market residual must be 3T-1"))
    core=Vector{T}(undef,2n-1)
    stage3_market_residual!(core,w,ctx;fixed=false)
    copyto!(view(res,1:2n-1),core)
    scale=max(one(T),inp.reference.Y)
    @inbounds for t in 1:n
        if t<n
            goods=w.Y[t]-w.I[t]-w.adjustment[t]-
                  inp.reference.export_shift*w.terms[t]^(-p.theta_star)
            household=w.aggregate_home[t]
        else
            last_gap=w.I[t-1]/w.K[t]-p.delta
            terminal_cost=inp.phi_ac/T(2)*w.K[t]*last_gap^2
            goods=w.Y[t]-w.I[t-1]-terminal_cost-
                  inp.reference.export_shift*w.terms[t]^(-p.theta_star)
            # Match period-T resources to period-T household spending.
            household=w.aggregate_home[t]
        end
        res[2n-1+t]=(goods-household)/scale
    end
    all(isfinite,res) ||
        throw(ArgumentError("nonfinite explicit flexible GE residual"))
    return res
end

function stage3_flexible_full_residual!(res::AbstractVector{T},
                                       x::AbstractVector{T},
                                       ctx::Stage3Context{T},
                                       w::TransitionWorkspace{T}) where {T<:AbstractFloat}
    w.function_evals+=1
    stage3_flexible_full_prices!(w,ctx,x)
    stage3_household_backward!(w,ctx)
    stage3_distribution_forward!(w,ctx)
    return stage3_flexible_full_market_residual!(res,w,ctx)
end

"""
    solve_stage3_flexible_full(ctx; x0, kwargs...)

Solve the 599-equation FULL flexible-price GE, without a nested scalar
trade-price root. x0 must contain a 399-element K/L guess followed by all
200 positive candidate terms of trade. Any returned result must still pass
an independent household-budget, probability and goods-market audit.
"""
function solve_stage3_flexible_full(ctx::Stage3Context{T};
                                    x0::AbstractVector{T},
                                    kwargs...) where {T<:AbstractFloat}
    return stage3_nonlinear_solve(ctx,Vector{T}(x0);
                                  explicit_flexible=true,kwargs...)
end
