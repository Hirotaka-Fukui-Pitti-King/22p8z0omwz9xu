"""
Firm-level state initialization and idiosyncratic shocks for Nirei (2015),
Sections 2.2, 2.6, 3.2–3.3, equation (16), and Appendix C.

There are no aggregate exogenous shocks. All firm-level draws are performed
with explicitly supplied RNG streams. The lambda profile is fixed in time.
"""

"""
    FirmRNGs

Separate pseudorandom streams for indivisibilities, initial productivity,
initial capital gaps and subsequent productivity innovations.
The stored UInt64 seeds are reproducibility metadata, not a mathematical
guarantee that different pseudorandom streams are independent.
"""
struct FirmRNGs
    lambda_seed::UInt64
    initial_productivity_seed::UInt64
    initial_gap_seed::UInt64
    subsequent_productivity_seed::UInt64
    lambda_rng::Xoshiro
    initial_productivity_rng::Xoshiro
    initial_gap_rng::Xoshiro
    subsequent_productivity_rng::Xoshiro
end

"""
    firm_rngs(base_seed=20261101, replication=1)

Construct four reproducible Xoshiro streams, varying deterministically
across Monte Carlo replications. Fix Julia and the Manifest to reproduce
the exact stream; Julia does not promise bit-for-bit RNG compatibility
across releases.
"""
function firm_rngs(base_seed::Integer=20261101, replication::Integer=1)
    base_seed >= 0 || throw(ArgumentError("base_seed must be nonnegative"))
    replication >= 1 || throw(ArgumentError("replication must be >= 1"))
    root = UInt64(base_seed) + UInt64(replication - 1) * 0x9e3779b97f4a7c15
    seed_lambda = root ⊻ 0x243f6a8885a308d3
    seed_a0 = root ⊻ 0x13198a2e03707344
    seed_gap = root ⊻ 0xa4093822299f31d0
    seed_at = root ⊻ 0x082efa98ec4e6c89
    return FirmRNGs(
        seed_lambda, seed_a0, seed_gap, seed_at,
        Xoshiro(seed_lambda), Xoshiro(seed_a0),
        Xoshiro(seed_gap), Xoshiro(seed_at),
    )
end

stream_seeds(rngs::FirmRNGs) = (
    lambda = rngs.lambda_seed,
    initial_productivity = rngs.initial_productivity_seed,
    initial_gap = rngs.initial_gap_seed,
    subsequent_productivity = rngs.subsequent_productivity_seed,
)

"""
    FirmState

All arrays have the same length. Current productivity is a_{i,t};
next_productivity is a preallocated buffer for a_{i,t+1}, initially
undefined until draw_productivity! is called. Subsequent investment
decisions use current capital and next productivity (Issue #4).
A mutable outer struct permits O(1) buffer swapping at period transitions.
The initially sampled gaps and thresholds are kept for documentation;
capital and productivity vectors are updated in place by later issues.
"""
mutable struct FirmState
    lambda::Vector{Float64}
    capital::Vector{Float64}
    productivity::Vector{Float64}
    next_productivity::Vector{Float64}
    initial_gap::Vector{Float64}
    initial_b::Vector{Float64}
    next_ready::Bool
end

"""
    threshold_coefficient(lambda, productivity, calibration)

Return b_i in equation (16), where k*_i=b_i*K. Use log/expm1 to
avoid cancellation for lambda near 1.
"""
function threshold_coefficient(
    lambda::Float64,
    productivity::Float64,
    calibration::CalibrationResult,
)
    lambda > 1.0 || throw(DomainError(lambda, "lambda must exceed one"))
    productivity > 0.0 && isfinite(productivity) ||
        throw(DomainError(productivity, "productivity must be positive and finite"))
    p = calibration.parameters
    rho = calibration.constants.rho
    A = calibration.threshold.A
    x = log(lambda)
    logratio = log(expm1(rho * x)) - log(expm1(x))
    return exp(log(A) + ((rho / p.alpha) * log(productivity) + logratio) /
                          (1.0 - rho))
end

"""
    aggregate_capital(productivity, capital, alpha, rho)

Finite-N productivity-weighted aggregate capital from equation (13).
The supplied productivity array must correspond to the same capital
vintage. No allocations proportional to N.
"""
function aggregate_capital(
    productivity::Vector{Float64},
    capital::Vector{Float64},
    alpha::Float64,
    rho::Float64,
)
    length(productivity) == length(capital) && !isempty(capital) ||
        throw(DimensionMismatch("firm arrays must have the same positive length"))
    0.0 < alpha < 1.0 && 0.0 < rho < 1.0 ||
        throw(ArgumentError("alpha and rho must lie in (0,1)"))
    weight = rho / alpha
    total = 0.0
    @inbounds for i in eachindex(productivity, capital)
        total += productivity[i]^weight * capital[i]^rho
    end
    K = (total / length(capital))^(inv(rho))
    isfinite(K) && K > 0.0 ||
        throw(DomainError(K, "finite-N aggregate capital must be positive and finite"))
    return K
end

function aggregate_capital(
    state::FirmState,
    calibration::CalibrationResult;
    use_next_productivity::Bool=false,
)
    use_next_productivity && !state.next_ready &&
        throw(ArgumentError("draw_productivity! must precede use_next_productivity=true"))
    a = use_next_productivity ? state.next_productivity : state.productivity
    return aggregate_capital(
        a, state.capital, calibration.parameters.alpha, calibration.constants.rho
    )
end

"""
    initialize_firms(calibration, rngs; N=calibration.parameters.firms)

Initialize N firms from Appendix C: lambda_i ~ truncated normal,
a_{i,0} independent lognormal, s_{i,0} ~ Uniform[0,1).
The finite sample is built with Kbar as reference:
k_{i,0} = A*[a_{i,0}^(rho/alpha)*(lambda_i^rho-1)/(lambda_i-1)]
                 ^(1/(1-rho))*Kbar*lambda_i^s_{i,0}.

This preserves each *sampled* gap relative to the reference Kbar, but
the actual finite-N K0 is calculated rather than forced to equal Kbar.
The paper specifies the distributions, not a finite-sample scaling rule.
"""
function initialize_firms(
    calibration::CalibrationResult,
    rngs::FirmRNGs;
    N::Integer=calibration.parameters.firms,
)
    N > 0 || throw(ArgumentError("N must be positive"))
    n = Int(N)
    p = calibration.parameters
    c = calibration.constants
    s = calibration.steady_state
    lambda = Vector{Float64}(undef, n)
    a0 = Vector{Float64}(undef, n)
    gap0 = Vector{Float64}(undef, n)
    b0 = Vector{Float64}(undef, n)
    k0 = Vector{Float64}(undef, n)

    # Distributions.jl scalar sampling; no temporary N-sized draw array.
    d = lambda_distribution(p)
    @inbounds for i in eachindex(lambda)
        lambda[i] = rand(rngs.lambda_rng, d)
    end
    randn!(rngs.initial_productivity_rng, a0)
    rand!(rngs.initial_gap_rng, gap0)

    @inbounds for i in eachindex(k0)
        ai = exp(c.mean_log_a + p.sigma_log_a * a0[i])
        a0[i] = ai
        b = threshold_coefficient(lambda[i], ai, calibration)
        b0[i] = b
        k0[i] = b * s.K * exp(gap0[i] * log(lambda[i]))
    end

    all(x -> isfinite(x) && x > 0.0, a0) ||
        throw(ArgumentError("invalid initial firm productivity"))
    all(x -> isfinite(x) && x > 0.0, k0) ||
        throw(ArgumentError("invalid initial firm capital"))
    return FirmState(
        lambda, k0, a0, Vector{Float64}(undef, n),
        gap0, b0, false,
    )
end

"""
    draw_productivity!(state, rngs, calibration)

Draw independent a_{i,t+1} into a preallocated buffer without modifying
lambda, current capital, or current productivity.
"""
function draw_productivity!(
    state::FirmState,
    rngs::FirmRNGs,
    calibration::CalibrationResult,
)
    a_next = state.next_productivity
    length(a_next) == length(state.productivity) ||
        throw(DimensionMismatch("next productivity buffer has wrong length"))
    c = calibration.constants
    sd = calibration.parameters.sigma_log_a
    randn!(rngs.subsequent_productivity_rng, a_next)
    @inbounds @simd for i in eachindex(a_next)
        a_next[i] = exp(c.mean_log_a + sd * a_next[i])
    end
    state.next_ready = true
    return a_next
end

"""
    advance_productivity!(state)

Promote the previously drawn a_{i,t+1} to the current state. Swap
array references in O(1); no new N-sized allocation. A later equilibrium
implementation must update capital to k_{i,t+1} before calling this.
"""
function advance_productivity!(state::FirmState)
    state.next_ready || throw(ArgumentError("next-period productivity not drawn"))
    previous = state.productivity
    state.productivity = state.next_productivity
    state.next_productivity = previous
    state.next_ready = false
    return state
end

"""
    realized_initial_gap_summary(state, calibration)

Report the gap implied by actual finite-N K0 rather than by the
reference Kbar. It need not lie in [0,1) if finite sampling shifts K0.
No projection, clipping, resampling, or hidden rescaling is performed.
"""
function realized_initial_gap_summary(
    state::FirmState,
    calibration::CalibrationResult,
)
    n = length(state.lambda)
    K0 = aggregate_capital(state, calibration)
    shift = log(calibration.steady_state.K / K0)
    s = 0.0
    s2 = 0.0
    lo = Inf
    hi = -Inf
    outside = 0
    @inbounds for i in eachindex(state.lambda)
        observed = state.initial_gap[i] + shift / log(state.lambda[i])
        s += observed
        s2 += observed^2
        lo = min(lo, observed)
        hi = max(hi, observed)
        outside += (observed < 0.0 || observed >= 1.0)
    end
    avg = s / n
    dispersion = sqrt(max(0.0, (s2 - n * avg^2) / max(n - 1, 1)))
    return (
        K0 = K0,
        Kbar = calibration.steady_state.K,
        K0_relative_deviation = K0 / calibration.steady_state.K - 1.0,
        realized_gap_mean = avg,
        realized_gap_std = dispersion,
        realized_gap_min = lo,
        realized_gap_max = hi,
        realized_gap_outside_0_1_count = outside,
    )
end
