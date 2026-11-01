"""
Exact finite-N, one-period investment equilibria for Nirei (2015),
equations (10), (13), (14), (16) and Section 3.4.

The reaction function is a nondecreasing STEP function. Sweep its 2N
breakpoints rather than iterating a fictitious tatonnement (ES2).
Original Section 3.4 selection minimizes |log(X_actual/X_expected)| over
positive-investment fixed points. Opt-in :capital_fallback extends Definition 1
only when Section 3.4's log objective has no real-domain candidate.
"""

"""
    EquilibriumWorkspace

Reusable O(N) buffers for scanning 2N events, holding all fixed-point
candidates, and computing the selected next-period capital vector.
The sweep does not modify the current FirmState.
"""
mutable struct EquilibriumWorkspace
    switches::Vector{Float64}       # odd = lower ell_i, even = upper u_i
    order::Vector{Int}              # sorted permutation of 1:2N
    depreciated::Vector{Float64}    # (1-delta)*k_i,t
    base_power::Vector{Float64}     # a_i,t+1^(rho/alpha)*d_i^rho/N
    next_capital::Vector{Float64}
    fixed_capital::Vector{Float64}
    fixed_investment::Vector{Float64}
    n_fixed::Int
end

function initialize_equilibrium_workspace(n::Integer)
    n > 0 || throw(ArgumentError("number of firms must be positive"))
    k = Int(n)
    return EquilibriumWorkspace(
        Vector{Float64}(undef, 2k), Vector{Int}(undef, 2k),
        Vector{Float64}(undef, k), Vector{Float64}(undef, k),
        Vector{Float64}(undef, k),
        Vector{Float64}(undef, 2k + 1), Vector{Float64}(undef, 2k + 1),
        0,
    )
end

initialize_equilibrium_workspace(state::FirmState) =
    initialize_equilibrium_workspace(length(state.capital))

"""
    FixedPointEnumeration

Summary of the exact interval sweep. Candidate capital and investment
arrays live in the reusable workspace, in their first n_fixed entries.
"""
struct FixedPointEnumeration
    total::Int
    admissible::Int
    nonpositive::Int
    minimum_K::Float64
    maximum_K::Float64
end

struct ES1Selection
    index::Int
    criterion::Float64
    expected_K::Float64
    admissible::Int
    nonpositive::Int
end

"""
    EquilibriumResult

Post-selection one-period result. On successful return, the FirmState
has already advanced from (k_t,a_t) to (k_{t+1},a_{t+1}) by O(1) buffer
swaps. K_next is the selected fixed point, while the residual records
its comparison against the exact N-firm aggregation.
"""
struct EquilibriumResult
    K_prior::Float64
    expected_K::Float64
    expected_X::Float64
    K_next::Float64
    X::Float64
    epsilon::Float64
    objective::Float64
    n_fixed::Int
    n_admissible::Int
    n_nonpositive::Int
    minimum_K::Float64
    maximum_K::Float64
    n_investing::Int
    n_inactive::Int
    n_divesting::Int
    capital_residual::Float64
    investment_residual::Float64
    max_threshold_violation::Float64
end

# Candidate K belongs to (left, right], as required by
# kappa=-1 at K=ell_i and kappa=0 at K=u_i. No tolerance is used to
# admit a candidate that lies outside the interval.
@inline function _record_interval!(
    workspace::EquilibriumWorkspace,
    power_sum::Float64,
    investment_sum::Float64,
    left::Float64,
    right::Float64,
    invrho::Float64,
)
    power_sum > 0.0 && isfinite(power_sum) ||
        throw(DomainError(power_sum, "aggregate power sum must be finite and positive"))
    candidate = power_sum^invrho
    isfinite(candidate) ||
        throw(DomainError(candidate, "candidate aggregate capital is nonfinite"))
    if left < candidate <= right
        j = workspace.n_fixed + 1
        workspace.fixed_capital[j] = candidate
        workspace.fixed_investment[j] = investment_sum
        workspace.n_fixed = j
    end
    return nothing
end

"""
    enumerate_fixed_points!(workspace, state, calibration)

Enumerate ALL finite-N fixed points K=Phi(K) from equation (13), using
a single O(N log N) sorting pass and O(N) cumulative scanning.

For K <= ell_i the firm divests (kappa=-1), for ell_i < K <= u_i
it is inactive (kappa=0), and for K > u_i it invests (kappa=+1).

Initially all firms divest. The scan tests interval (left,right]
BEFORE applying events at right. Events with exactly equal breakpoint
values are applied simultaneously AFTER the interval is checked.
Floating point tolerances are not used in interval membership.
The current firm capital/productivity arrays are not mutated.
"""
function enumerate_fixed_points!(
    workspace::EquilibriumWorkspace,
    state::FirmState,
    calibration::CalibrationResult,
)
    state.next_ready ||
        throw(ArgumentError("call draw_productivity! before equilibrium enumeration"))
    p = calibration.parameters
    rho = calibration.constants.rho
    n = length(state.lambda)
    n > 0 && length(state.capital) == n &&
        length(state.next_productivity) == n ||
        throw(DimensionMismatch("FirmState arrays differ in length"))
    length(workspace.depreciated) == n ||
        throw(DimensionMismatch("EquilibriumWorkspace size must equal firm count"))

    invn = inv(Float64(n))
    invrho = inv(rho)
    power_sum = 0.0
    investment_sum = 0.0

    @inbounds for i in 1:n
        lambda = state.lambda[i]
        a = state.next_productivity[i]
        d = (1.0 - p.delta) * state.capital[i]
        lambda > 1.0 && isfinite(lambda) && isfinite(a) && a > 0.0 &&
            isfinite(d) && d > 0.0 ||
            throw(ArgumentError("invalid firm state at index $i"))

        b = threshold_coefficient(lambda, a, calibration)
        lower = d / (lambda * b)
        upper = d / b
        isfinite(lower) && isfinite(upper) && 0.0 < lower < upper ||
            throw(DomainError((lower, upper), "invalid switching thresholds for firm $i"))
        base = (a^(rho / p.alpha)) * d^rho * invn
        isfinite(base) && base > 0.0 ||
            throw(DomainError(base, "invalid firm power contribution $i"))

        workspace.depreciated[i] = d
        workspace.base_power[i] = base
        workspace.switches[2i - 1] = lower
        workspace.switches[2i] = upper
        power_sum += base * lambda^(-rho)
        investment_sum += d * (inv(lambda) - 1.0) * invn
    end

    sortperm!(workspace.order, workspace.switches)
    workspace.n_fixed = 0
    left = 0.0
    j = 1
    m = 2n
    @inbounds while j <= m
        boundary = workspace.switches[workspace.order[j]]
        _record_interval!(workspace, power_sum, investment_sum,
                          left, boundary, invrho)

        # At K = boundary, the *old* policies still apply. After
        # crossing it, process the complete batch before the next interval.
        while j <= m && workspace.switches[workspace.order[j]] == boundary
            event = workspace.order[j]
            firm = (event + 1) >>> 1
            lambda = state.lambda[firm]
            d = workspace.depreciated[firm]
            base = workspace.base_power[firm]
            if isodd(event)      # ell_i: -1 -> 0
                power_sum += base * (1.0 - lambda^(-rho))
                investment_sum += d * (1.0 - inv(lambda)) * invn
            else                 # u_i: 0 -> +1
                power_sum += base * (lambda^rho - 1.0)
                investment_sum += d * (lambda - 1.0) * invn
            end
            j += 1
        end
        left = boundary
    end
    _record_interval!(workspace, power_sum, investment_sum,
                      left, Inf, invrho)

    total = workspace.n_fixed
    total > 0 || throw(ErrorException("No finite-firm fixed point was found"))
    admissible = 0
    minK = Inf
    maxK = -Inf
    @inbounds for i in 1:total
        K = workspace.fixed_capital[i]
        minK = min(minK, K)
        maxK = max(maxK, K)
        admissible += (workspace.fixed_investment[i] > 0.0)
    end
    return FixedPointEnumeration(total, admissible,
                                 total - admissible, minK, maxK)
end

"""
    select_es1(workspace, expected_X, expected_K)

Choose the minimum |log X - log expected_X| over ALL enumerated fixed
points having X>0 (Section 3.4). Nonpositive X is counted and excluded
because its real logarithm is undefined; no clipping is performed.
Opt-in selection_mode=:capital_fallback preserves the original choice
whenever any X>0 root exists, but selects min |log(K/expected_K)| among
ALL roots when none do (Section 2.3, Definition 1). This fallback is NOT
specified by Section 3.4 and must be reported as a separate experiment.
When original scores tie exactly, select closest capital and then smallest K.
"""
function select_es1(
    workspace::EquilibriumWorkspace,
    expected_X::Real,
    expected_K::Real;
    selection_mode::Symbol=:original,
)
    selection_mode in (:original, :capital_fallback) ||
        throw(ArgumentError("selection_mode must be :original or :capital_fallback"))
    xe, ke = Float64(expected_X), Float64(expected_K)
    isfinite(xe) && xe > 0.0 ||
        throw(DomainError(xe, "expected investment must be positive and finite"))
    isfinite(ke) && ke > 0.0 ||
        throw(DomainError(ke, "expected capital must be positive and finite"))
    workspace.n_fixed > 0 ||
        throw(ArgumentError("enumerate_fixed_points! must run before select_es1"))

    best = 0
    best_score = Inf
    best_distance = Inf
    positive = 0
    @inbounds for i in 1:workspace.n_fixed
        X = workspace.fixed_investment[i]
        if X > 0.0
            positive += 1
            score = abs(log(X / xe))
            distance = abs(log(workspace.fixed_capital[i] / ke))
            if score < best_score ||
                (score == best_score &&
                 (distance < best_distance ||
                  (distance == best_distance &&
                   (best == 0 || workspace.fixed_capital[i] < workspace.fixed_capital[best]))))
                best = i
                best_score = score
                best_distance = distance
            end
        end
    end
    if positive == 0 && selection_mode == :capital_fallback
        # Definition 1 capital-distance completion only when Section 3.4
        # has no positive-investment candidates. Physical signed X is kept.
        @inbounds for i in 1:workspace.n_fixed
            K = workspace.fixed_capital[i]
            score = abs(log(K / ke))
            if score < best_score ||
                (score == best_score &&
                 (best == 0 || K < workspace.fixed_capital[best]))
                best = i
                best_score = score
            end
        end
    elseif positive == 0
        # Diagnostic on the genuine finite-N equilibrium correspondence.
        # Do not clip X or silently switch to capital-based selection:
        # Section 3.4 defines the ES1 shock using log(X/X_expected).
        fixed_ks = view(workspace.fixed_capital, 1:workspace.n_fixed)
        fixed_xs = view(workspace.fixed_investment, 1:workspace.n_fixed)
        throw(DomainError((
            fixed_points=workspace.n_fixed,
            minimum_investment=minimum(fixed_xs),
            maximum_investment=maximum(fixed_xs),
            minimum_capital=minimum(fixed_ks),
            maximum_capital=maximum(fixed_ks),
            expected_investment=xe,
            expected_capital=ke,
        ), "No equilibrium has X>0, so the Section 3.4 ES1 log shock is undefined; " *
           "the capital-based ES1 of Definition 1 is not a valid substitute for log(X/X_expected)"))
    end
    return ES1Selection(best, best_score, ke,
                        positive, workspace.n_fixed - positive)
end

"""
    solve_equilibrium!(workspace, state, calibration, expected_X, K_prior)

Solve one exact finite-N period after draw_productivity! has produced
a_{i,t+1}. Expected_X comes from Issue #2; K_prior must be the SAME
predetermined aggregate capital used in forming that expectation.

Steps: enumerate all fixed points, select ES1, form next firm capital
in workspace, independently check aggregation and net investment,
then atomically swap next capital and productivity into FirmState.

The expected capital reference for tie-breaking or the opt-in negative-X
fallback is continuum equation (18): (1-delta)*K_prior + (A^(rho-1)/rho)*expected_X.
This is distinct from the finite-N realization of aggregate capital.

Only the selected fixed point is reaggregated over N firms, never all
candidate fixed points. The original mode is unchanged. In the opt-in
fallback, a selected X<=0 retains its signed level; its epsilon is a
FIRST-ORDER capital-innovation proxy log(K/Kstar)/delta, NOT the undefined
Section 3.4 log-investment surprise.
"""
function solve_equilibrium!(
    workspace::EquilibriumWorkspace,
    state::FirmState,
    calibration::CalibrationResult,
    expected_X::Real,
    K_prior::Real;
    selection_mode::Symbol=:original,
)
    p = calibration.parameters
    rho = calibration.constants.rho
    kp = Float64(K_prior)
    xe = Float64(expected_X)
    isfinite(kp) && kp > 0.0 ||
        throw(DomainError(kp, "K_prior must be positive and finite"))
    isfinite(xe) && xe > 0.0 ||
        throw(DomainError(xe, "expected_X must be positive and finite"))

    fixed_points = enumerate_fixed_points!(workspace, state, calibration)
    Kstar = (1.0 - p.delta) * kp +
            calibration.threshold.A^(rho - 1.0) / rho * xe
    selection = select_es1(workspace, xe, Kstar; selection_mode=selection_mode)
    kselected = workspace.fixed_capital[selection.index]
    xselected = workspace.fixed_investment[selection.index]
    n = length(state.capital)

    investing = 0
    inactive = 0
    divesting = 0
    aggregate_investment = 0.0
    worst_gap = 0.0

    @inbounds for i in 1:n
        ell = workspace.switches[2i - 1]
        upper = workspace.switches[2i]
        lambda = state.lambda[i]
        d = workspace.depreciated[i]
        if kselected <= ell
            newk = d / lambda
            divesting += 1
        elseif kselected <= upper
            newk = d
            inactive += 1
        else
            newk = d * lambda
            investing += 1
        end
        workspace.next_capital[i] = newk
        aggregate_investment += (newk - d)
        worst_gap = max(worst_gap,
                        kselected <= ell ? max(0.0, kselected - ell) :
                        kselected <= upper ? max(0.0, ell - kselected, kselected - upper) :
                        max(0.0, upper - kselected))
    end

    K_check = aggregate_capital(
        state.next_productivity, workspace.next_capital, p.alpha, rho
    )
    X_check = aggregate_investment / n
    K_residual = K_check - kselected
    X_residual = X_check - xselected
    K_scale = max(1.0, abs(kselected))
    X_scale = max(1.0, abs(xselected))
    abs(K_residual) <= 1e-10 * K_scale ||
        throw(ErrorException("finite-N K fixed-point residual = $K_residual"))
    abs(X_residual) <= 1e-10 * X_scale ||
        throw(ErrorException("finite-N investment residual = $X_residual"))
    investing + inactive + divesting == n ||
        throw(ErrorException("firm policy counts do not sum to N"))
    (selection_mode == :capital_fallback || xselected > 0.0) ||
        throw(DomainError(xselected, "selected original Section 3.4 ES1 investment must be positive"))
    # When X<=0, log(X/Xe) does not exist. The capital forecast
    # innovation / delta is a separately labelled first-order surrogate.
    innovation = xselected > 0.0 ? log(xselected / xe) :
                 log(kselected / Kstar) / p.delta

    result = EquilibriumResult(
        kp, Kstar, xe, kselected, xselected, innovation,
        selection.criterion, fixed_points.total, selection.admissible,
        selection.nonpositive, fixed_points.minimum_K, fixed_points.maximum_K,
        investing, inactive, divesting, K_residual, X_residual, worst_gap,
    )
    # Commit the next period's state only after all numerical checks pass.
    previous_capital = state.capital
    state.capital = workspace.next_capital
    workspace.next_capital = previous_capital
    advance_productivity!(state)
    return result
end
