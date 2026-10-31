"""
    BenchmarkParameters

Quarterly benchmark from Nirei (2015), Table 1 and Appendix C.
The standard deviation of lambda is that of the *untruncated* normal.
Firm-level draws are deferred to Issue #3.
"""
Base.@kwdef struct BenchmarkParameters
    alpha::Float64 = 0.26
    delta::Float64 = 0.02
    beta::Float64 = 0.99
    sigma::Float64 = 15.0
    eta::Float64 = 11.0
    psi::Float64 = 1.0
    zeta::Float64 = 2.0
    lambda_mean::Float64 = 1.028
    lambda_sd::Float64 = 0.004
    lambda_truncation_sd::Float64 = 2.0
    sigma_log_a::Float64 = 0.0005
    firms::Int = 110_000
end

struct StructuralConstants
    c_L::Float64
    rho::Float64
    productivity_exponent::Float64
    mean_log_a::Float64
end

struct ThresholdConstant
    A::Float64
    integral::Float64
    integral_error::Float64
end

struct SteadyState
    K::Float64
    Y::Float64
    C::Float64
    X::Float64
    L::Float64
    w::Float64
    R::Float64
    mu::Float64
    K_over_Y::Float64
    a_c::Float64
    a_x::Float64
    a_R::Float64
end

struct CalibrationResult
    parameters::BenchmarkParameters
    constants::StructuralConstants
    threshold::ThresholdConstant
    steady_state::SteadyState
end

benchmark_parameters() = BenchmarkParameters()

function _check_parameters(p::BenchmarkParameters)
    for (label, value) in (
        ("alpha", p.alpha), ("delta", p.delta), ("beta", p.beta),
        ("sigma", p.sigma), ("eta", p.eta), ("psi", p.psi),
        ("zeta", p.zeta), ("lambda_mean", p.lambda_mean),
        ("lambda_sd", p.lambda_sd), ("lambda_truncation_sd", p.lambda_truncation_sd),
        ("sigma_log_a", p.sigma_log_a),
    )
        isfinite(value) || throw(ArgumentError("$label must be finite"))
    end
    0.0 < p.alpha < 1.0 || throw(ArgumentError("alpha must be in (0,1)"))
    0.0 < p.delta < 1.0 || throw(ArgumentError("delta must be in (0,1)"))
    0.0 < p.beta < 1.0 || throw(ArgumentError("beta must be in (0,1)"))
    p.sigma > 0.0 || throw(ArgumentError("sigma must be positive"))
    p.eta > 1.0 || throw(ArgumentError("eta must exceed one"))
    p.psi > 0.0 || throw(ArgumentError("psi must be positive"))
    p.zeta > 1.0 || throw(ArgumentError("zeta must exceed one"))
    p.lambda_mean > 1.0 || throw(ArgumentError("lambda_mean must exceed one"))
    p.lambda_sd > 0.0 || throw(ArgumentError("lambda_sd must be positive"))
    p.lambda_truncation_sd > 0.0 ||
        throw(ArgumentError("lambda_truncation_sd must be positive"))
    p.sigma_log_a > 0.0 || throw(ArgumentError("sigma_log_a must be positive"))
    p.firms > 0 || throw(ArgumentError("firms must be positive"))
    lower = p.lambda_mean - p.lambda_truncation_sd * p.lambda_sd
    lower > 1.0 || throw(ArgumentError("lambda distribution requires lower bound > 1"))
    return nothing
end

"""
    lambda_distribution(p)

Normal(lambda_mean, lambda_sd) truncated at ±lambda_truncation_sd original
standard deviations. This is a distribution only; no firm draws are made.
"""
function lambda_distribution(p::BenchmarkParameters)
    _check_parameters(p)
    lower = p.lambda_mean - p.lambda_truncation_sd * p.lambda_sd
    upper = p.lambda_mean + p.lambda_truncation_sd * p.lambda_sd
    return truncated(Normal(p.lambda_mean, p.lambda_sd), lower, upper)
end

"""
    compute_structural_constants(p)

Compute c_L, rho, q=rho/[alpha*(1-rho)], and the lognormal mean that
normalizes E[a^q]=1. All calculations are deterministic.
"""
function compute_structural_constants(p::BenchmarkParameters)
    _check_parameters(p)
    c_L = (1.0 - inv(p.eta)) * (1.0 - p.alpha)
    rho = (1.0 - inv(p.eta)) * p.alpha / (1.0 - c_L)
    0.0 < rho < 1.0 || throw(DomainError(rho, "rho must lie in (0,1)"))
    q = rho / (p.alpha * (1.0 - rho))
    mean_log_a = -0.5 * q * p.sigma_log_a^2
    return StructuralConstants(c_L, rho, q, mean_log_a)
end

# Numerically stable evaluation of the lambda integrand in equation (16).
# Writing the first factor in log form avoids fractional powers of small
# quantities and cancellation in lambda^rho - 1.
function _threshold_integrand(lambda::Float64, rho::Float64)
    x = log(lambda)
    logfactor = (log(expm1(rho * x)) - rho * log(expm1(x))) / (1.0 - rho)
    return exp(logfactor) / (rho * x)
end

"""
    compute_threshold_constant(p, constants; rtol=1e-11)

Numerically compute A in Nirei (2015), equation (16). Appendix C sets
E[a^(rho/(alpha*(1-rho)))]=1; hence the remaining expectation is over
the truncated-normal lambda distribution only. Return the quadrature
estimate and QuadGK's absolute integration error estimate.
"""
function compute_threshold_constant(
    p::BenchmarkParameters,
    constants::StructuralConstants;
    rtol::Float64 = 1e-11,
)
    0.0 < rtol < 1.0 || throw(ArgumentError("rtol must be in (0,1)"))
    distribution = lambda_distribution(p)
    lower = p.lambda_mean - p.lambda_truncation_sd * p.lambda_sd
    upper = p.lambda_mean + p.lambda_truncation_sd * p.lambda_sd
    rho = constants.rho

    integrand(lambda) = _threshold_integrand(lambda, rho) * pdf(distribution, lambda)
    integral, error = quadgk(integrand, lower, upper; rtol=rtol)

    isfinite(integral) && integral > 0.0 ||
        throw(DomainError(integral, "equation (16) expectation must be positive and finite"))
    isfinite(error) && error >= 0.0 ||
        throw(DomainError(error, "quadrature error estimate must be finite"))
    A = integral^(-inv(rho))
    isfinite(A) && A > 0.0 || throw(DomainError(A, "A must be positive and finite"))
    return ThresholdConstant(A, integral, error)
end

"""
    compute_steady_state(p, constants, threshold)

Compute the continuum steady state of equations (7)–(12), (17), and (18).
No general-equilibrium fixed-point iteration is needed at this stage.
"""
function compute_steady_state(
    p::BenchmarkParameters,
    constants::StructuralConstants,
    threshold::ThresholdConstant,
)
    alpha, delta, beta = p.alpha, p.delta, p.beta
    c_L, rho, A = constants.c_L, constants.rho, threshold.A

    R = inv(beta)
    rental_denom = R - 1.0 + delta
    a_x = delta * rho * (1.0 - c_L) / rental_denom
    a_c = 1.0 - a_x
    0.0 < a_x < 1.0 ||
        throw(DomainError(a_x, "steady-state investment share must be in (0,1)"))

    K_over_Y = (1.0 - c_L) / (A^(1.0 - rho) * rental_denom)
    L = (c_L / (p.psi * (c_L + p.zeta * a_c)))^(inv(p.zeta))
    # Production Y = K^alpha * L^(1-alpha), K = (K/Y)*Y.
    Y = K_over_Y^(alpha / (1.0 - alpha)) * L
    K = K_over_Y * Y
    X = a_x * Y
    C = a_c * Y
    w = c_L * Y / L
    leisure_term = 1.0 - p.psi * L^p.zeta
    leisure_term > 0.0 ||
        throw(DomainError(leisure_term, "utility requires 1 - psi*L^zeta > 0"))
    mu = C^(-p.sigma) * leisure_term^(1.0 - p.sigma)
    a_R = R / rental_denom

    vals = (K, Y, C, X, L, w, R, mu, K_over_Y, a_c, a_x, a_R)
    all(isfinite, vals) || throw(ArgumentError("steady state contains nonfinite values"))
    all(>(0.0), vals) || throw(ArgumentError("steady state requires positive values"))
    return SteadyState(vals...)
end

"""
    equilibrium_residuals(result)

Minimal embedded consistency checks of the steady-state equations.
Returns relative or normalized residuals; no separate diagnostic project.
"""
function equilibrium_residuals(result::CalibrationResult)
    p = result.parameters
    c = result.constants
    t = result.threshold
    s = result.steady_state
    return (
        resource_market = (s.C + s.X - s.Y) / s.Y,
        production = (s.K^p.alpha * s.L^(1.0 - p.alpha) - s.Y) / s.Y,
        labor_share = (s.w * s.L / s.Y - c.c_L),
        labor_supply = (s.w * (1.0 - p.psi * s.L^p.zeta) -
                        s.C * p.psi * p.zeta * s.L^(p.zeta - 1.0)) / s.w,
        capital_accumulation = (p.delta * s.K - t.A^(c.rho - 1.0) * s.X / c.rho) / s.K,
        euler = p.beta * s.R - 1.0,
        productivity_moment = exp(c.productivity_exponent * c.mean_log_a +
            0.5 * (c.productivity_exponent * p.sigma_log_a)^2) - 1.0,
    )
end

"""
    run_calibration(; parameters=benchmark_parameters(), rtol=1e-11)

Return all components used later by Issues #2 and #3. Abort if the
original equilibrium identities fail to floating-point accuracy.
"""
function run_calibration(;
    parameters::BenchmarkParameters=benchmark_parameters(),
    rtol::Float64=1e-11,
)
    c = compute_structural_constants(parameters)
    t = compute_threshold_constant(parameters, c; rtol=rtol)
    s = compute_steady_state(parameters, c, t)
    result = CalibrationResult(parameters, c, t, s)
    residuals = equilibrium_residuals(result)
    max_residual = maximum(abs, Tuple(residuals))
    max_residual <= 1e-10 ||
        throw(ErrorException("steady-state equation residual = $max_residual"))
    return result
end
