"""
One-trial dynamic simulation for Nirei (2015), Sections 3.3-3.5 and
Appendix A. All firm-level decisions and the finite-N equilibrium are
provided by the already implemented Issues #3 and #4.

Timeline: before period t's investment decision, current firm capital
and productivity aggregate to K_t; mu_{t-1} is predetermined. Appendix A
forecasts (K_{t+1}, mu_t, X_t, ...) from (K_t, mu_{t-1}). New a_{i,t+1}
are drawn, ES1 realizes X_t and K_{t+1}, then macro variables for date t
are reconstructed. The new state has (K_{t+1}, mu_t).

FIRST-ORDER CONTINUATION (not a solved nonlinear finite-N RE equilibrium):
Appendix A forecasts log(mu_t) conditionally and gives the capital innovation
delta*epsilon_t. The local stable eigenvector slope h maps that innovation
to the current marginal-utility surprise. The realized exact finite-N K
may differ from the first-order predicted K; this discrepancy is recorded.
Realized Y,C,L,w and R satisfy exact equations (7),(8),(11),(12) and the
marginal-utility definition, on the locally continuous household branch.
This fixes erroneous exponentiation of the linearized resource equation.
Neither step asserts global nonlinear rational-expectations consistency.
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
Realized Y,C,L,w,R are exact level allocations conditional on the first-order
marginal-utility continuation. For backward-compatible CSV field names, the
five *_linear_residual columns now store NORMALIZED EXACT LEVEL residuals,
not Appendix A log-linear residuals. Their values are zero at a valid root.
The mismatch between realized finite-N capital and Appendix A linear capital
is recorded separately.
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

# Conditional on the (first-order) marginal utility realization, solve the
# EXACT contemporaneous household/production/market-clearing equations.
# Use the low-labor branch, which contains the published deterministic steady
# state; the second (high-labor) root is not silently selected. Neither the
# predetermined firm-price expectation equations (9),(17) nor a realized
# constant labor-share FOC are imposed as false ex-post identities.
#
# For X>0 and sigma>1 the log marginal-utility function in L has ONE
# minimum: L*derivative(log mu)= -sigma*(1-alpha)*Y/C
#   + (sigma-1)*psi*zeta*L^zeta/(1-psi*L^zeta).
# Both terms are monotone increasing in L, and have opposite infinite
# endpoint signs. Thus the low-L root is unique whenever it exists.
function _realized_macroeconomy(
    calibration::CalibrationResult,
    K_prior::Float64,
    X::Float64,
    logmu::Float64,
    logmu_prior::Float64,
)
    p = calibration.parameters
    s = calibration.steady_state
    isfinite(K_prior) && K_prior > 0.0 ||
        throw(DomainError(K_prior, "K_prior must be positive"))
    isfinite(X) && X > 0.0 ||
        throw(DomainError(X, "original Section 3.4 ES1 requires X>0"))
    p.sigma > 1.0 ||
        throw(DomainError(p.sigma, "low-labor branch proof requires sigma>1"))
    isfinite(logmu) && isfinite(logmu_prior) ||
        throw(DomainError((logmu, logmu_prior), "nonfinite marginal utility state"))

    capital_factor = K_prior^p.alpha
    labor_exponent = 1.0 - p.alpha
    Lmax = p.psi^(-inv(p.zeta))
    Lmin = (X / capital_factor)^(inv(labor_exponent))
    isfinite(Lmin) && Lmin < Lmax ||
        throw(DomainError((K_prior, X), "no labor with positive C and feasible leisure"))
    span = Lmax - Lmin
    lo = Lmin + 1e-12 * span
    hi = Lmax - 1e-12 * span
    lo < hi || throw(DomainError((Lmin, Lmax), "empty numerical labor bracket"))

    function logmu_labor(L::Float64)
        Y = capital_factor * L^labor_exponent
        C = Y - X
        leisure = 1.0 - p.psi * L^p.zeta
        C > 0.0 && leisure > 0.0 ||
            throw(DomainError((L, C, leisure), "household allocation outside domain"))
        return -p.sigma * log(C) + (1.0 - p.sigma) * log(leisure)
    end
    function scaled_derivative(L::Float64)
        Y = capital_factor * L^labor_exponent
        C = Y - X
        leisure = 1.0 - p.psi * L^p.zeta
        return -p.sigma * labor_exponent * Y / C +
               (p.sigma - 1.0) * p.psi * p.zeta * L^p.zeta / leisure
    end

    # Locate the minimum of log(mu(L)) on the feasible labor interval.
    # This determines whether the steady-state-connected low-labor branch
    # exists, but is NOT itself a stable way to invert extremely large mu:
    # at large mu, C=Y-X loses relative precision as C approaches zero.
    a, b = lo, hi
    for _ in 1:100
        mid = a + (b - a) / 2.0
        if scaled_derivative(mid) < 0.0
            a = mid
        else
            b = mid
        end
    end
    Lturn = a + (b - a) / 2.0
    target_logmu = log(s.mu) + logmu
    minimum_logmu = logmu_labor(Lturn)
    target_logmu >= minimum_logmu - 1e-12 ||
        throw(DomainError((target_logmu, minimum_logmu),
                          "no exact level household allocation for prescribed mu,K,X"))

    # Invert marginal utility algebraically instead of obtaining extremely
    # small consumption as the difference of two nearly equal goods flows.
    # mu = C^(-sigma) * (1-psi*L^zeta)^(1-sigma), hence
    # log(C)=[(1-sigma)*log(1-psi*L^zeta)-log(mu)]/sigma.
    function consumption_given_labor(L::Float64)
        leisure = 1.0 - p.psi * L^p.zeta
        leisure > 0.0 ||
            throw(DomainError((L, leisure), "nonpositive leisure"))
        C_implied = exp(((1.0 - p.sigma) * log(leisure) - target_logmu) / p.sigma)
        C_implied > 0.0 && isfinite(C_implied) ||
            throw(DomainError(C_implied, "consumption unresolvable at Float64 precision"))
        return C_implied
    end
    function goods_gap(L::Float64)
        return capital_factor * L^labor_exponent - X -
               consumption_given_labor(L)
    end

    # The lower branch has negative goods_gap at Lmin (zero goods
    # available for C), and positive goods_gap at the mu minimum
    # whenever the target mu is strictly above that minimum.
    L = Lturn
    if target_logmu > minimum_logmu + 1e-12
        a, b = Lmin, Lturn
        goods_gap(a) <= 0.0 ||
            throw(ErrorException("cannot bracket low-labor root at Float64 precision"))
        goods_gap(b) >= -1e-12 * max(1.0, abs(X)) ||
            throw(ErrorException("cannot bracket low-labor root near mu minimum"))
        for _ in 1:120
            mid = a + (b - a) / 2.0
            if mid == a || mid == b
                break
            elseif goods_gap(mid) < 0.0
                a = mid
            else
                b = mid
            end
        end
        # Choose the representable endpoint with the smaller resource
        # residual, rather than taking another midpoint that may round.
        L = abs(goods_gap(a)) <= abs(goods_gap(b)) ? a : b
    end
    Y = capital_factor * L^labor_exponent
    C = consumption_given_labor(L)
    leisure = 1.0 - p.psi * L^p.zeta
    w = C * p.psi * p.zeta * L^(p.zeta - 1.0) / leisure
    mu = s.mu * exp(logmu)
    mu_prior = s.mu * exp(logmu_prior)
    R = mu_prior / (p.beta * mu)

    all(z -> isfinite(z) && z > 0.0, (Y, C, L, w, R, mu, mu_prior)) ||
        throw(DomainError((Y, C, L, w, R), "nonpositive exact macro level"))
    production_res = (Y - K_prior^p.alpha * L^labor_exponent) / Y
    resources_res = (Y - C - X) / Y
    mu_res = (logmu_labor(L) - target_logmu)
    labor_res = (w * leisure - C * p.psi * p.zeta * L^(p.zeta - 1.0)) / w
    euler_res = p.beta * R * mu / mu_prior - 1.0
    maximum(abs, (production_res, resources_res, mu_res,
                  labor_res, euler_res)) < 1e-9 ||
        throw(ErrorException("exact level household/market equations failed"))

    return (
        Y=Y, C=C, L=L, w=w, R=R,
        logY=log(Y / s.Y), logC=log(C / s.C),
        logL=log(L / s.L), logw=log(w / s.w),
        logR=log(R / s.R),
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

The marginal utility surprise follows the LOCAL first-order Appendix A
stable response to the log-investment innovation: logmu_t = Eprev logmu_t
+ h*delta*epsilon_t, rather than mechanically imposing h*REALIZED finite-N
log capital. These differ when Eq(26)'s approximation for K is inaccurate.
This is NOT a fully solved global nonlinear rational expectations system.
The realized household/market-clearing levels are solved exactly and are
rejected if no admissible low-labor branch exists.
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
        # Appendix A: the capital innovation in the LINEAR system is
        # delta*epsilon, not the discrepancy in exact finite-N CES capital.
        # The stable eigenspace slope maps this innovation to the
        # marginal-utility jump (first-order RE closure).
        logmu = expected.logmu_next + state.re.h * p.delta * eq.epsilon
        mu = s.mu * exp(logmu)
        isfinite(mu) && mu > 0.0 ||
            throw(DomainError(mu, "unbounded/nonpositive continuation marginal utility"))
        realized_macro = _realized_macroeconomy(
            state.calibration, kp, eq.X, logmu, logmu_prior
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
            realized_macro.Y, realized_macro.C, eq.X, realized_macro.L, realized_macro.w, realized_macro.R, mu,
            log(eq.K_next), log(realized_macro.Y), log(realized_macro.C),
            log(eq.X), log(realized_macro.L),
            logK_prior, logK_next, logmu,
            expected.K_next, expected.X_next, expected.C_next,
            expected.L_next, expected.Y_next, expected.mu_next,
            expected.w_next, expected.R_next,
            eq.epsilon, logK_error, logmu_error, stable_jump_residual,
            eq.n_fixed, eq.n_admissible, eq.n_nonpositive, eq.objective,
            eq.n_investing, eq.n_inactive, eq.n_divesting,
            eq.capital_residual, eq.investment_residual,
            realized_macro.production_residual, realized_macro.resource_residual,
            realized_macro.mu_residual, realized_macro.labor_residual,
            realized_macro.euler_residual, (realized_macro.Y - realized_macro.C - eq.X) / realized_macro.Y,
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
