"""
One-trial dynamic simulation for Nirei (2015), Sections 3.3-3.5 and
Appendix A. All firm-level decisions and the finite-N equilibrium are
provided by the already implemented Issues #3 and #4.

Timeline: before period t's investment decision, current firm capital
and productivity aggregate to K_t; mu_{t-1} is predetermined. Appendix A
forecasts (K_{t+1}, mu_t, X_t, ...) from (K_t, mu_{t-1}). New a_{i,t+1}
are drawn, ES1 realizes X_t and K_{t+1}, then macro variables for date t
are reconstructed. The new state has (K_{t+1}, mu_t).

CLOSURE (not an independently estimated finite-N rational-expectations
solution): imposing the bounded/stable continuation on the forward-looking
marginal utility gives tilde(mu)_t = h*tilde(K)_{t+1}.
This is the saddle-path terminal/transversality restriction, not an
equation of predetermined prices at realized outcomes. It agrees with
the conditional RE system *to first order* when investment innovations
have zero conditional mean. Finite-N ES1 need not have that property:
the resulting forecast errors and their drift are recorded, NOT removed.
"""

"""
    SimulationState

Only one period of firm arrays and the single reusable 2N-event workspace
are stored. The linear RE eigensystem is constructed just once.
"""
mutable struct SimulationState
    calibration::CalibrationResult
    re::LinearRESolution
    firms::FirmState
    rngs::FirmRNGs
    workspace::EquilibriumWorkspace
    K::Float64
    mu_prior::Float64
    period::Int
    base_seed::Int
    replication::Int
end

struct SimulationFailure <: Exception
    period::Int
    K_prior::Float64
    mu_prior::Float64
    reason::String
end

function Base.showerror(io::IO, ex::SimulationFailure)
    print(io, "Nirei2015 simulation failed at t=", ex.period,
          " (K_prior=", ex.K_prior, ", mu_prior=", ex.mu_prior, "): ", ex.reason)
end

"""
    PeriodRecord

Log columns are log *levels* (suitable for later Table 2), while
logdev_* columns are deviations from the fixed calibration steady state.
All realized C,L,Y,w,mu,R levels are exp-transforms of the FIRST-ORDER
macro system, not solutions of its nonlinear levels counterpart.
"""
struct PeriodRecord
    t::Int
    K_prior::Float64
    K_next_firm::Float64
    K_next_linear::Float64
    capital_linear_gap::Float64
    Y::Float64
    C::Float64
    X::Float64
    L::Float64
    w::Float64
    R::Float64
    mu::Float64
    logK::Float64
    logY::Float64
    logC::Float64
    logX::Float64
    logL::Float64
    logdev_K_prior::Float64
    logdev_K_next::Float64
    logdev_mu::Float64
    expected_K_next::Float64
    expected_X::Float64
    expected_C::Float64
    expected_L::Float64
    expected_Y::Float64
    expected_mu::Float64
    expected_w::Float64
    expected_R::Float64
    epsilon::Float64
    logK_expectation_error::Float64
    logmu_expectation_error::Float64
    stable_jump_residual::Float64
    n_fixed::Int
    n_positive_X::Int
    n_nonpositive_X::Int
    es1_objective::Float64
    n_investing::Int
    n_inactive::Int
    n_divesting::Int
    fixed_point_residual::Float64
    investment_residual::Float64
    production_linear_residual::Float64
    resources_linear_residual::Float64
    marginal_utility_linear_residual::Float64
    labor_supply_linear_residual::Float64
    euler_linear_residual::Float64
    level_resources_gap::Float64
end

"""
    initialize_simulation(; calibration=run_calibration(), base_seed=20261101,
                            replication=1, N=calibration.parameters.firms)

Initialize firms using Issue #3. Initialize mu_{-1} on the bounded
stable eigenspace. Subsequently the period transition uses the same
bounded-continuation condition to determine the forward-looking mu_t.
The finite-N initial aggregate capital is never rescaled to Kbar.
"""
function initialize_simulation(;
    calibration::CalibrationResult=run_calibration(),
    base_seed::Int=20261101,
    replication::Int=1,
    N::Int=calibration.parameters.firms,
)
    re = compute_linear_re_solution(calibration)
    rngs = firm_rngs(base_seed, replication)
    firms = initialize_firms(calibration, rngs; N=N)
    workspace = initialize_equilibrium_workspace(firms)
    K = aggregate_capital(firms, calibration)
    mu_prior = calibration.steady_state.mu *
               exp(re.h * log(K / calibration.steady_state.K))
    isfinite(mu_prior) && mu_prior > 0.0 ||
        throw(DomainError(mu_prior, "nonpositive initial marginal utility"))
    return SimulationState(calibration, re, firms, rngs, workspace,
                           K, mu_prior, 0, base_seed, replication)
end

# Return log deviations of realized C,L,Y,w and their linear-equation
# residuals. eqs (27),(28),(30),(33) are contemporaneous conditions.
# Price-setting conditions (29),(31) are NOT imposed on realizations.
function _realized_macroeconomy(
    calibration::CalibrationResult,
    re::LinearRESolution,
    logK_prior::Float64,
    logX::Float64,
    logmu::Float64,
    logmu_prior::Float64,
)
    p = calibration.parameters
    s = calibration.steady_state
    c = calibration.constants
    A1 = re.coefficients.A1
    abs(A1) > 1e-12 || throw(DomainError(A1, "singular macro closure"))
    logL = (s.a_x * logX - p.alpha * logK_prior -
            (s.a_c / p.sigma) * logmu) / A1
    logY = p.alpha * logK_prior + (1.0 - p.alpha) * logL
    logC = (logY - s.a_x * logX) / s.a_c
    logw = logC + (p.zeta - 1.0 + c.c_L / s.a_c) * logL
    # The realized inverse stochastic-discount-factor identity (7)/(32).
    # This is not an ex-ante interest-rate forecasting equation.
    logR = logmu_prior - logmu

    production_res = logY - p.alpha * logK_prior -
                     (1.0 - p.alpha) * logL
    resources_res = logY - s.a_c * logC - s.a_x * logX
    mu_res = logmu + p.sigma * logC -
             (p.sigma - 1.0) * (c.c_L / s.a_c) * logL
    labor_res = logw - logC -
                (p.zeta - 1.0 + c.c_L / s.a_c) * logL
    euler_res = logR - logmu_prior + logmu
    maxres = maximum(abs, (production_res, resources_res, mu_res,
                           labor_res, euler_res))
    maxres <= 1e-9 ||
        throw(ErrorException("realized linearized macro conditions failed: $maxres"))

    vals = (s.Y * exp(logY), s.C * exp(logC), s.L * exp(logL),
            s.w * exp(logw), s.R * exp(logR))
    all(x -> isfinite(x) && x > 0.0, vals) ||
        throw(ArgumentError("macro levels must be positive and finite"))
    p.psi * vals[3]^p.zeta < 1.0 ||
        throw(DomainError(vals[3], "level reconstructed labor violates utility domain"))
    return (
        Y=vals[1], C=vals[2], L=vals[3], w=vals[4], R=vals[5],
        logY=logY, logC=logC, logL=logL, logw=logw, logR=logR,
        production_residual=production_res,
        resource_residual=resources_res,
        mu_residual=mu_res,
        labor_residual=labor_res,
        euler_residual=euler_res,
    )
end

"""
    simulate_period!(state)

Perform one period in the order
  conditional expectations -> independent productivity draws ->
  exact firm equilibrium with ES1 -> realized macro variables ->
  continuation-state update.

The value mu_t is the contemporaneously adjusting forward-looking
marginal utility. We use the stable continuation relation h*K_{t+1}
as an explicit first-order selection assumption. It is NOT obtained
by treating predetermined expected firm prices as realized identities.
The departure of logmu_t from the actual conditional RE forecast is
retained as an endogenous forecast error. Replacing this continuation
condition with an independently derived nonlinear solution would
require a stronger solution concept not specified in Appendix A.
"""
function simulate_period!(state::SimulationState)
    p = state.calibration.parameters
    s = state.calibration.steady_state
    t = state.period
    kp = state.K
    mup = state.mu_prior
    try
        expected = compute_expected_equilibrium(
            state.re, state.calibration, kp, mup
        )
        draw_productivity!(state.firms, state.rngs, state.calibration)
        eq = solve_equilibrium!(
            state.workspace, state.firms, state.calibration, expected.X_next, kp
        )
        # Firm arrays already represent (k_{t+1},a_{t+1}) here.
        logK_prior = log(kp / s.K)
        logK_next = log(eq.K_next / s.K)
        logmu_prior = log(mup / s.mu)
        logmu = state.re.h * logK_next
        mu = s.mu * exp(logmu)
        logX = log(eq.X / s.X)
        isfinite(mu) && mu > 0.0 ||
            throw(DomainError(mu, "unbounded/nonpositive continuation marginal utility"))
        macro = _realized_macroeconomy(
            state.calibration, state.re, logK_prior, logX, logmu, logmu_prior
        )

        # Paper (26) is a FIRST-ORDER aggregate capital equation; it
        # must not replace the exact finite-firm equation (13).
        K_linear_logdev = expected.logK_next + p.delta * eq.epsilon
        K_linear = s.K * exp(K_linear_logdev)
        isfinite(K_linear) && K_linear > 0.0 ||
            throw(DomainError(K_linear, "invalid log-linear implied capital"))
        logK_error = logK_next - expected.logK_next
        logmu_error = logmu - expected.logmu_next
        stable_jump_residual = logmu_error - state.re.h * logK_error
        record = PeriodRecord(
            t, kp, eq.K_next, K_linear, logK_next - K_linear_logdev,
            macro.Y, macro.C, eq.X, macro.L, macro.w, macro.R, mu,
            log(eq.K_next), log(macro.Y), log(macro.C),
            log(eq.X), log(macro.L),
            logK_prior, logK_next, logmu,
            expected.K_next, expected.X_next, expected.C_next,
            expected.L_next, expected.Y_next, expected.mu_next,
            expected.w_next, expected.R_next,
            eq.epsilon, logK_error, logmu_error, stable_jump_residual,
            eq.n_fixed, eq.n_admissible, eq.n_nonpositive, eq.objective,
            eq.n_investing, eq.n_inactive, eq.n_divesting,
            eq.capital_residual, eq.investment_residual,
            macro.production_residual, macro.resource_residual,
            macro.mu_residual, macro.labor_residual,
            macro.euler_residual, (macro.Y - macro.C - eq.X) / macro.Y,
        )
        # The same finite-N capital drives the next period's expectations.
        state.K = eq.K_next
        state.mu_prior = mu
        state.period += 1
        return record
    catch err
        err isa InterruptException && rethrow()
        throw(SimulationFailure(t, kp, mup, sprint(showerror, err)))
    end
end

"""
    simulate_path!(state, T=1100)

Return exactly T consecutive PeriodRecord values. No full firm-panel
history is stored. On an error, throw SimulationFailure with the date and
predetermined state; the invoking CLI writes a failure summary and stops.
"""
function simulate_path!(state::SimulationState, T::Integer=1100)
    T > 0 || throw(ArgumentError("T must be positive"))
    records = Vector{PeriodRecord}(undef, Int(T))
    for j in eachindex(records)
        records[j] = simulate_period!(state)
    end
    return records
end
