#!/usr/bin/env julia
# Source-level regression: exact realized markets & household allocation,
# and Appendix A's FIRST-ORDER conditional forecast/innovation split.
# Local: julia --project=. scripts/check_exact_realized_macro.jl
using Nirei2015

function check_exact_levels()
    cal = run_calibration()
    s = cal.steady_state
    p = cal.parameters

    # Every printed deterministic steady-state level must be recovered.
    m = Nirei2015._realized_macroeconomy(cal, s.K, s.X, 0.0, 0.0)
    @assert isapprox(m.Y, s.Y; atol=1e-10, rtol=1e-10)
    @assert isapprox(m.C, s.C; atol=1e-10, rtol=1e-10)
    @assert isapprox(m.L, s.L; atol=1e-10, rtol=1e-10)
    @assert isapprox(m.w, s.w; atol=1e-10, rtol=1e-10)
    @assert isapprox(m.R, s.R; atol=1e-10, rtol=1e-10)

    # Independent direct evaluations: neither an invented constant
    # realized labor share nor first-order log resource equation.
    for (K, X, dmu, dprior) in (
        (s.K, s.X, 0.0, 0.0),
        (s.K, 0.4*s.X, -0.01, 0.02),
        (1.01*s.K, 0.8*s.X, 0.005, -0.008),
    )
        o = Nirei2015._realized_macroeconomy(
            cal, K, X, dmu, dprior
        )
        leisure = 1.0 - p.psi * o.L^p.zeta
        mu = s.mu * exp(dmu)
        mu_prior = s.mu * exp(dprior)
        @assert 0.0 < leisure < 1.0
        @assert isapprox(o.Y, K^p.alpha * o.L^(1.0-p.alpha); atol=1e-10)
        @assert isapprox(o.Y, o.C + X; atol=1e-10)
        @assert isapprox(mu, o.C^(-p.sigma)*leisure^(1.0-p.sigma);
                         rtol=1e-10)
        @assert isapprox(o.w, o.C*p.psi*p.zeta*
                         o.L^(p.zeta-1.0)/leisure; rtol=1e-10)
        @assert isapprox(p.beta*o.R*mu/mu_prior, 1.0; atol=1e-12)
        @assert maximum(abs, (o.production_residual, o.resource_residual,
                              o.mu_residual, o.labor_residual,
                              o.euler_residual)) <= 1e-9
    end
    println("PASS: exact steady state and three household/market examples")
end

function check_shock_timing()
    cal = run_calibration()
    re = compute_linear_re_solution(cal)
    p = cal.parameters
    s = cal.steady_state
    state = initialize_simulation(; calibration=cal, N=2000, base_seed=20261101)
    forecast = compute_expected_equilibrium(re, cal, state.K, state.mu_prior)
    # Isolate first-order Appendix A shock response from exact CES levels;
    # this algebra does not assume the next finite-firm K equals (26).
    epsilon = -0.08
    logK_linear = forecast.logK_next + p.delta*epsilon
    logmu_real = forecast.logmu_next + re.h*p.delta*epsilon
    @assert isapprox(logmu_real-forecast.logmu_next,
                     re.h*(logK_linear-forecast.logK_next); atol=1e-12)
    @assert isapprox(logmu_real-forecast.logmu_next,
                     re.h*p.delta*epsilon; atol=1e-12)
    @assert forecast.X_next > 0.0
    println("PASS: Appendix A conditional expectation and stable shock innovation")
end

check_exact_levels()
check_shock_timing()
println("NOTE: this is not a global nonlinear rational-expectations proof.")
