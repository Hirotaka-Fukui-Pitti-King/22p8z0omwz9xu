"""
    LinearRECoefficients

Coefficients derived from Nirei (2015), Appendix A, equations (26)-(37).
"""
struct LinearRECoefficients
    eta_L::Float64
    A1::Float64
    a_R::Float64
end

"""
    LinearRESolution

Conditional-mean transition z_next = F*z_prior, where

    z_prior = [tilde(K)_{-1}, tilde(mu)_{-1}]
    z_next  = [E_{-1} tilde(K)_0, E_{-1} tilde(mu)_0].

The first component of the next-period realized capital differs from the
forecast by delta*epsilon_0 (Appendix A, after equation (37)).
All eigenvalue and matrix calculations are done once, not inside a time loop.
"""
struct LinearRESolution
    coefficients::LinearRECoefficients
    B::Matrix{Float64}
    D::Matrix{Float64}
    F::Matrix{Float64}
    eigenvalues::Vector{Float64}
    stable_index::Int
    stable_eigenvalue::Float64
    unstable_eigenvalue::Float64
    stable_eigenvector::Vector{Float64}
    h::Float64
    matrix_residual::Float64
    eigen_residual::Float64
    saddle_residual::Float64
end

"""
    ExpectedEquilibrium

All tilde quantities are conditional *expected log deviations* at date 0,
conditional on the information observable at date -1. Fields ending in
_next represent t=0, and fields ending in _prior represent t=-1.
The reported levels are exp-transformed first-order log forecasts,
not exact conditional means of levels.
"""
struct ExpectedEquilibrium
    K_prior::Float64
    mu_prior::Float64
    logK_prior::Float64
    logmu_prior::Float64
    logK_next::Float64
    logmu_next::Float64
    logX_next::Float64
    logC_next::Float64
    logL_next::Float64
    logY_next::Float64
    logw_next::Float64
    logR_next::Float64
    K_next::Float64
    mu_next::Float64
    X_next::Float64
    C_next::Float64
    L_next::Float64
    Y_next::Float64
    w_next::Float64
    R_next::Float64
end

"""
    compute_linear_re_solution(calibration)

Obtain the unique saddle-path solution of Nirei (2015), Proposition 5,
under the benchmark condition a_x < alpha. B and D are reproduced from
equations (36) and (37), with B21=D12=0. Reject non-real eigenvalues and
configurations without exactly one stable and one unstable eigenvalue.
"""
function compute_linear_re_solution(calibration::CalibrationResult)
    p = calibration.parameters
    c = calibration.constants
    s = calibration.steady_state

    s.a_x < p.alpha ||
        throw(DomainError(s.a_x, "strict saddle stability requires a_x < alpha"))
    eta_L = inv((2.0 - inv(p.sigma)) * (c.c_L / s.a_c) + p.zeta - 1.0)
    A1 = 1.0 - p.alpha - c.c_L * (p.sigma - 1.0) / p.sigma
    a_R = s.a_R
    isfinite(eta_L) && eta_L > 0.0 ||
        throw(DomainError(eta_L, "Frisch elasticity must be positive"))
    coefficients = LinearRECoefficients(eta_L, A1, a_R)

    # Equation (37): B11*E(K0) + B12*E(mu0) = D11*K_-1.
    B11 = s.a_x / p.delta
    B12 = -s.a_c / p.sigma -
          A1 / (p.sigma * (p.alpha + inv(eta_L)))
    D11 = p.alpha + s.a_x * (1.0 - p.delta) / p.delta +
          A1 * p.alpha / (p.alpha + inv(eta_L))

    # Equation (36): B22*E(mu0) = D21*K_-1 + D22*mu_-1.
    B22 = eta_L / (p.sigma * (1.0 + p.alpha * eta_L)) +
          a_R / (1.0 - p.alpha)
    D21 = inv(1.0 + p.alpha * eta_L)
    D22 = a_R / (1.0 - p.alpha)

    B = [B11 B12; 0.0 B22]
    D = [D11 0.0; D21 D22]
    F = B \ D
    eigen_decomposition = eigen(F)
    complex_values = eigen_decomposition.values
    maximum(abs, imag.(complex_values)) <= 1e-10 ||
        throw(DomainError(complex_values, "expected system has non-real eigenvalues"))
    eigenvalues = Float64[real(x) for x in complex_values]

    stable = findall(x -> abs(x) < 1.0 - 1e-9, eigenvalues)
    unstable = findall(x -> abs(x) > 1.0 + 1e-9, eigenvalues)
    length(stable) == 1 && length(unstable) == 1 ||
        throw(DomainError(eigenvalues, "expected exactly one stable and one unstable root"))
    stable_index = only(stable)
    xi_s = eigenvalues[stable_index]
    xi_u = eigenvalues[only(unstable)]

    raw_v = eigen_decomposition.vectors[:, stable_index]
    maximum(abs, imag.(raw_v)) <= 1e-10 ||
        throw(DomainError(raw_v, "stable eigenvector must be real"))
    abs(real(raw_v[1])) > 1e-12 ||
        throw(DomainError(raw_v, "stable eigenvector capital component too small"))
    v = Float64[real(raw_v[1]), real(raw_v[2])]
    v ./= v[1]
    h = v[2]

    matrix_residual = norm(B * F - D, Inf)
    eigen_residual = norm(F * v - xi_s * v, Inf)
    saddle_residual = abs(F[2, 1] + F[2, 2] * h - h * xi_s)
    maximum((matrix_residual, eigen_residual, saddle_residual)) < 1e-10 ||
        throw(ErrorException("Appendix A linear system residual exceeds 1e-10"))

    return LinearRESolution(
        coefficients, B, D, F, eigenvalues, stable_index,
        xi_s, xi_u, v, h,
        matrix_residual, eigen_residual, saddle_residual,
    )
end

"""
    compute_expected_equilibrium(solution, calibration, K_prior, mu_prior)

One-step conditional expected equilibrium in the Appendix A information timing:
K_prior=K_{-1} and mu_prior=mu_{-1} are both predetermined/observable
before the period-0 investment shock. mu_prior is NOT silently reset
to h*log(K_prior/Kbar), since a shock may have displaced the state.
"""
function compute_expected_equilibrium(
    solution::LinearRESolution,
    calibration::CalibrationResult,
    K_prior::Real,
    mu_prior::Real,
)
    s = calibration.steady_state
    p = calibration.parameters
    c = calibration.constants
    kp = Float64(K_prior)
    mup = Float64(mu_prior)
    isfinite(kp) && kp > 0.0 || throw(DomainError(kp, "K_prior must be positive and finite"))
    isfinite(mup) && mup > 0.0 ||
        throw(DomainError(mup, "mu_prior must be positive and finite"))

    logkp = log(kp / s.K)
    logmup = log(mup / s.mu)
    F = solution.F
    logkn = F[1, 1] * logkp + F[1, 2] * logmup
    logmun = F[2, 1] * logkp + F[2, 2] * logmup

    # (26), (35), (33), (27), (30), (32).
    logxn = (logkn - (1.0 - p.delta) * logkp) / p.delta
    logln = (p.alpha * logkp + logmun / p.sigma) /
            (p.alpha + inv(solution.coefficients.eta_L))
    logcn = (-logmun + (p.sigma - 1.0) * (c.c_L / s.a_c) * logln) / p.sigma
    logyn = p.alpha * logkp + (1.0 - p.alpha) * logln
    logwn = logcn + (p.zeta - 1.0 + c.c_L / s.a_c) * logln
    logrn = logmup - logmun

    levels = (
        s.K * exp(logkn), s.mu * exp(logmun), s.X * exp(logxn),
        s.C * exp(logcn), s.L * exp(logln),
        s.Y * exp(logyn), s.w * exp(logwn), s.R * exp(logrn),
    )
    all(isfinite, levels) && all(x -> x > 0.0, levels) ||
        throw(ArgumentError("expected equilibrium contains nonfinite/nonpositive levels"))

    return ExpectedEquilibrium(
        kp, mup, logkp, logmup, logkn, logmun,
        logxn, logcn, logln, logyn, logwn, logrn,
        levels...,
    )
end

"""
    compute_saddle_path_expectations(solution, calibration, K_prior)

Convenience for a *no-new-shock stable-path state*, where
tilde(mu)_{-1} = h*tilde(K)_{-1}. Not a permissible automatic
normalization of the realized state after a stochastic shock.
"""
function compute_saddle_path_expectations(
    solution::LinearRESolution,
    calibration::CalibrationResult,
    K_prior::Real,
)
    kp = Float64(K_prior)
    s = calibration.steady_state
    isfinite(kp) && kp > 0.0 || throw(DomainError(kp, "K_prior must be positive and finite"))
    mu_prior = s.mu * exp(solution.h * log(kp / s.K))
    return compute_expected_equilibrium(solution, calibration, kp, mu_prior)
end

"""
    expectation_residuals(solution, calibration, expected)

Residuals of each of equations (26)-(33), (35)-(37), written at t=-1/0.
Expected variables are conditional expectations of log-deviations.
"""
function expectation_residuals(
    solution::LinearRESolution,
    calibration::CalibrationResult,
    expected::ExpectedEquilibrium,
)
    p = calibration.parameters
    s = calibration.steady_state
    c = calibration.constants
    co = solution.coefficients
    k, mu = expected.logK_prior, expected.logmu_prior
    kn, mun = expected.logK_next, expected.logmu_next
    x, cc, l = expected.logX_next, expected.logC_next, expected.logL_next
    y, w, r = expected.logY_next, expected.logw_next, expected.logR_next

    bstate = solution.B * [kn, mun]
    dstate = solution.D * [k, mu]

    return (
        system = norm(bstate - dstate, Inf),
        eq26_capital = kn - ((1.0 - p.delta) * k + p.delta * x),
        eq27_output = y - (p.alpha * k + (1.0 - p.alpha) * l),
        eq28_resources = y - (s.a_c * cc + s.a_x * x),
        eq29_firm_wage = y - (k - (1.0 - p.alpha) / p.alpha * w),
        eq30_labor_supply = w - (cc + (p.zeta - 1.0 + c.c_L / s.a_c) * l),
        eq31_factor_price = (1.0 - p.alpha) / p.alpha * w + s.a_R * r,
        eq32_euler = r - (mu - mun),
        eq33_marginal_utility = mun - (-p.sigma * cc +
                                           (p.sigma - 1.0) * (c.c_L / s.a_c) * l),
        eq35_labor = -mun / p.sigma - p.alpha * k +
                     (p.alpha + inv(co.eta_L)) * l,
        eq36_mu = co.eta_L / (p.sigma * (1.0 + p.alpha * co.eta_L)) * mun +
                  co.a_R / (1.0 - p.alpha) * mun -
                  k / (1.0 + p.alpha * co.eta_L) -
                  co.a_R / (1.0 - p.alpha) * mu,
        eq37_capital = s.a_x / p.delta * kn -
                       (s.a_c / p.sigma +
                        co.A1 / (p.sigma * (p.alpha + inv(co.eta_L)))) * mun -
                       (p.alpha + s.a_x * (1.0 - p.delta) / p.delta +
                        co.A1 * p.alpha / (p.alpha + inv(co.eta_L))) * k,
    )
end
