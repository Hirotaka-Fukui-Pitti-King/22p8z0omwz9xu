#!/usr/bin/env julia
# Finite-N mechanism audit for the first investment collapse; not a model alternative.
# From Nirei-2015/: julia --project=. scripts/audit_investment_dynamics.jl
using Nirei2015
using TOML
using Distributions: Normal, cdf, ccdf

function response_at_K(firms, calibration, K::Float64; use_prior_productivity::Bool=false)
    p = calibration.parameters
    rho = calibration.constants.rho
    n = length(firms.capital)
    power = 0.0
    investment = 0.0
    ninvest = 0
    ninactive = 0
    ndivest = 0
    @inbounds for i in 1:n
        d = (1.0 - p.delta) * firms.capital[i]
        l = firms.lambda[i]
        a = use_prior_productivity ? firms.productivity[i] : firms.next_productivity[i]
        b = threshold_coefficient(l, a, calibration)
        nextk = if K <= d / (l*b)
            ndivest += 1
            d / l
        elseif K <= d / b
            ninactive += 1
            d
        else
            ninvest += 1
            d*l
        end
        power += a^(rho/p.alpha) * nextk^rho
        investment += nextk - d
    end
    return (
        Phi=(power/n)^(1/rho),
        X=investment/n,
        investing=ninvest,
        inactive=ninactive,
        divesting=ndivest,
    )
end

function gap_diagnostics(state)
    firms, calibration = state.firms, state.calibration
    n = length(firms.capital)
    total = 0.0
    total2 = 0.0
    outsider = 0
    high = -Inf
    low = Inf
    # Current capital and productivity use the same time index.
    @inbounds for i in 1:n
        b = threshold_coefficient(
            firms.lambda[i], firms.productivity[i], calibration
        )
        g = log(firms.capital[i]/(b*state.K)) / log(firms.lambda[i])
        low = min(low, g)
        high = max(high, g)
        outsider += (g < 0.0 || g >= 1.0)
        total += g
        total2 += g^2
    end
    m = total/n
    return Dict(
        "mean_gap"=>m,
        "sd_gap"=>sqrt(max(0.0, total2/n-m^2)),
        "min_gap"=>low, "max_gap"=>high,
        "fraction_gap_outside_0_1"=>outsider/n,
    )
end


# In the model's one-step policy, the gap after adjustment need not remain
# in [0,1) when the preferred log threshold moves by more than one
# indivisibility interval. Proposition 4 uses a modulo-1 evolution and
# the paper's footnote explicitly discusses expanding the choice set
# when productivity support is too wide. This is a DIAGNOSTIC only:
# no modification to the firm's choice set or to the benchmark shock law.
function threshold_support_diagnostics(state, K::Float64)
    firms, calibration = state.firms, state.calibration
    n = length(firms.capital)
    p = calibration.parameters
    cutoff_above_one = 0
    cutoff_below_zero = 0
    next_gap_outside = 0
    next_gap_below_zero = 0
    next_gap_above_one = 0
    change_sum = 0.0
    change_squared = 0.0
    for i in 1:n
        l = firms.lambda[i]
        b0 = threshold_coefficient(l, firms.productivity[i], calibration)
        b1 = threshold_coefficient(l, firms.next_productivity[i], calibration)
        ln_l = log(l)
        s0 = log(firms.capital[i] / (b0 * state.K)) / ln_l
        q = -log1p(-p.delta) / ln_l
        log_b_change = log(b1 / b0)
        cutoff = q + (log_b_change + log(K/state.K)) / ln_l
        cutoff_above_one += cutoff > 1.0
        cutoff_below_zero += cutoff < 0.0
        # Retain exactly the current {-1,0,1} policy.
        adjustment = s0 < cutoff ? 1 : s0 >= cutoff + 1.0 ? -1 : 0
        next_gap = s0 - cutoff + adjustment
        next_gap_outside += next_gap < 0.0 || next_gap >= 1.0
        next_gap_below_zero += next_gap < 0.0
        next_gap_above_one += next_gap >= 1.0
        change_sum += log_b_change
        change_squared += log_b_change^2
    end
    average = change_sum / n
    return Dict(
        "fraction_cutoff_above_one" => cutoff_above_one/n,
        "fraction_cutoff_below_zero" => cutoff_below_zero/n,
        "fraction_one_step_next_gap_outside_0_1" => next_gap_outside/n,
        "fraction_one_step_next_gap_below_zero" => next_gap_below_zero/n,
        "fraction_one_step_next_gap_above_one" => next_gap_above_one/n,
        "sd_log_threshold_change" =>
            sqrt(max(0.0, change_squared/n-average^2)),
        "mean_log_threshold_change" => average,
    )
end


# Explicit counterfactual based on Nirei (2015) Section 2.1 footnote 2
# and Section 2.6 footnote 4: allow the discrete adjustment exponent to
# take as many integer steps as needed to place the new gap in [0,1).
# THIS IS NOT the three-action benchmark; it changes the feasible action
# set and must never silently replace the published {-1,0,+1} policy.
function multi_step_response_at_K(firms, calibration, K::Float64)
    p = calibration.parameters
    rho = calibration.constants.rho
    n = length(firms.capital)
    power = 0.0
    investment = 0.0
    nmore_positive = 0
    nmore_negative = 0
    kmin = typemax(Int)
    kmax = typemin(Int)
    max_gap_violation = 0.0
    @inbounds for i in 1:n
        d = (1.0 - p.delta) * firms.capital[i]
        l = firms.lambda[i]
        a = firms.next_productivity[i]
        b = threshold_coefficient(l, a, calibration)
        # At exact equality, ceiling implements the original convention:
        # K=d/(lambda*b) -> -1; K=d/b -> 0; K=d*lambda/b -> 1.
        exponent = ceil(Int, (log(b*K) - log(d)) / log(l))
        exponent > 1 && (nmore_positive += 1)
        exponent < -1 && (nmore_negative += 1)
        kmin = min(kmin, exponent)
        kmax = max(kmax, exponent)
        nextk = d * l^exponent
        gap = (log(nextk) - log(b*K)) / log(l)
        max_gap_violation = max(
            max_gap_violation, max(0.0, -gap, gap - 1.0)
        )
        power += a^(rho/p.alpha) * nextk^rho
        investment += nextk - d
    end
    return (
        Phi=(power/n)^(1/rho),
        X=investment/n,
        n_more_than_one_investment=nmore_positive,
        n_more_than_one_divestment=nmore_negative,
        min_exponent=kmin,
        max_exponent=kmax,
        max_gap_violation=max_gap_violation,
    )
end


# Assumption audit: does not alter productivity, capital, the investment
# policy, the ES1 choice or the finite-N equilibrium enumeration.
#
# Section 2.1, p.957: |log(k*_{t+1}/[(1-delta)k*_t])| <= log(lambda_i).
# Section 2.6, p.964: Pr(|Delta log(a)| < q*alpha*(1/rho-1)*
#                        log(lambda_lower)) = 1.
# q in that passage references the anticipated investment fraction
# introduced in Section 2.4. Its analytical argument assumes q<0.5.
# Therefore the q=0.5 threshold below is an *upper-envelope diagnostic*,
# not a claim that the GE benchmark has q=0.5. The stationary q proxy is
# separately reported and is NOT substituted into the paper's assumption.
# Section 3.2, p.968: lambda_i*(1-delta)>1.
#
# The current benchmark draws an unbounded Normal for log(a), so no positive
# finite bound can hold with probability one. Per-period counts and the
# contribution of out-of-one-step-range firms remain empirical questions.
function boundedness_assumption_audit(
    state, K::Float64, expected_X::Float64, at_expected,
    multistep, stationary_q_proxy::Float64,
)
    firms, calibration = state.firms, state.calibration
    p = calibration.parameters
    rho = calibration.constants.rho
    n = length(firms.capital)
    ln_depreciation = log1p(-p.delta)
    lower_lambda = p.lambda_mean - p.lambda_truncation_sd*p.lambda_sd
    q_upper = 0.5
    paper_upper_bound = q_upper*p.alpha*(inv(rho)-1.0)*log(lower_lambda)
    stationary_proxy_bound =
        stationary_q_proxy*p.alpha*(inv(rho)-1.0)*log(lower_lambda)
    std_delta_log_a = sqrt(2.0)*p.sigma_log_a
    delta_normal = Normal(0.0, std_delta_log_a)
    analytic_p_upper = 2.0*ccdf(delta_normal,paper_upper_bound)
    analytic_p_proxy = 2.0*ccdf(delta_normal,stationary_proxy_bound)

    n_violate_s21 = 0
    n_violate_s26_upper = 0
    n_violate_s26_proxy = 0
    n_violate_s26_firm_specific = 0
    n_violate_s32 = 0
    n_multi = 0
    n_overlap_s21_multi = 0
    n_overlap_s26_multi = 0
    n_overlap_s32_multi = 0
    n_overlap_s21_s26 = 0
    n_pos_multiple = 0
    n_neg_multiple = 0
    logb_min = Inf
    logb_max = -Inf
    logb_sum = 0.0
    logb_sum_sq = 0.0
    threshold_movement_min = Inf
    threshold_movement_max = -Inf
    total_power_3 = 0.0
    extra_power = 0.0
    extra_investment = 0.0
    extra_power_on_s21 = 0.0
    extra_investment_on_s21 = 0.0
    extra_investment_on_s26 = 0.0
    extra_investment_on_s32 = 0.0

    @inbounds for i in 1:n
        lambda_i = firms.lambda[i]
        a_old = firms.productivity[i]
        a_new = firms.next_productivity[i]
        d = (1.0-p.delta)*firms.capital[i]
        b_old = threshold_coefficient(lambda_i,a_old,calibration)
        b_new = threshold_coefficient(lambda_i,a_new,calibration)
        ln_lambda = log(lambda_i)
        ln_a_change = log(a_new/a_old)
        ln_b_change = (rho/(p.alpha*(1.0-rho)))*ln_a_change
        threshold_move = ln_b_change + log(K/state.K)-ln_depreciation

        violation_s21 = abs(threshold_move)>ln_lambda
        violation_s26 = abs(ln_a_change)>=paper_upper_bound
        violation_s26_proxy = abs(ln_a_change)>=stationary_proxy_bound
        violation_s26_i = abs(ln_a_change)>=
            q_upper*p.alpha*(inv(rho)-1.0)*ln_lambda
        violation_s32 = lambda_i*(1.0-p.delta)<=1.0

        one_exp = if K<=d/(lambda_i*b_new)
            -1
        elseif K<=d/b_new
            0
        else
            1
        end
        ideal_exp = ceil(Int,log(b_new*K/d)/ln_lambda)
        extra = ideal_exp != one_exp
        n_violate_s21 += violation_s21
        n_violate_s26_upper += violation_s26
        n_violate_s26_proxy += violation_s26_proxy
        n_violate_s26_firm_specific += violation_s26_i
        n_violate_s32 += violation_s32
        n_multi += extra
        n_overlap_s21_multi += violation_s21 && extra
        n_overlap_s26_multi += violation_s26 && extra
        n_overlap_s32_multi += violation_s32 && extra
        n_overlap_s21_s26 += violation_s21 && violation_s26
        n_pos_multiple += ideal_exp > 1
        n_neg_multiple += ideal_exp < -1

        k3 = d*lambda_i^one_exp
        kideal = d*lambda_i^ideal_exp
        delta_x = (kideal-k3)/n
        delta_power_i = a_new^(rho/p.alpha)*(kideal^rho-k3^rho)/n
        total_power_3 += a_new^(rho/p.alpha)*k3^rho/n
        extra_power += delta_power_i
        extra_investment += delta_x
        if violation_s21
            extra_power_on_s21 += delta_power_i
            extra_investment_on_s21 += delta_x
        end
        if violation_s26
            extra_investment_on_s26 += delta_x
        end
        if violation_s32
            extra_investment_on_s32 += delta_x
        end
        logb_min = min(logb_min,ln_b_change)
        logb_max = max(logb_max,ln_b_change)
        threshold_movement_min = min(threshold_movement_min,threshold_move)
        threshold_movement_max = max(threshold_movement_max,threshold_move)
        logb_sum += ln_b_change
        logb_sum_sq += ln_b_change^2
    end

    mean_logb = logb_sum/n
    phi3 = total_power_3^(inv(rho))
    phi_with_extra = (total_power_3+extra_power)^(inv(rho))
    # Independent CES accounting identity: all increments are calculated
    # at fixed expected capital, never along a different equilibrium path.
    isapprox(phi3,at_expected.Phi;atol=1e-9,rtol=1e-10) ||
        error("original-policy CES decomposition disagrees with reaction")
    isapprox(phi_with_extra,multistep.Phi;atol=1e-9,rtol=1e-10) ||
        error("hypothetical CES decomposition disagrees with reaction")
    isapprox(extra_investment,multistep.X-at_expected.X;
             atol=1e-9,rtol=1e-10) ||
        error("investment decomposition disagrees with reaction")

    output = Dict{String,Any}(
        "q_analytical_upper_envelope" => q_upper,
        "q_stationary_proxy_not_paper_q" => stationary_q_proxy,
        "paper_s26_bound_using_q_upper" => paper_upper_bound,
        "s26_bound_using_stationary_proxy_for_comparison" =>
            stationary_proxy_bound,
        "s26_normal_difference_sd" => std_delta_log_a,
        "s26_analytical_violation_probability_q_upper" => analytic_p_upper,
        "s26_analytical_expected_count_q_upper" => n*analytic_p_upper,
        "s26_analytical_violation_probability_stationary_proxy" => analytic_p_proxy,
        "s26_analytical_expected_count_stationary_proxy" => n*analytic_p_proxy,
        "s21_count_expected_capital" => n_violate_s21,
        "s26_count_global_q_upper" => n_violate_s26_upper,
        "s26_count_global_stationary_proxy" => n_violate_s26_proxy,
        "s26_count_firm_specific_q_upper_not_paper_bound" =>
            n_violate_s26_firm_specific,
        "s32_count_lambda_depreciation" => n_violate_s32,
        "one_step_vs_unrestricted_different_action_count" => n_multi,
        "n_ideal_exponent_gt1" => n_pos_multiple,
        "n_ideal_exponent_lt_minus1" => n_neg_multiple,
        "s21_and_action_mismatch_count" => n_overlap_s21_multi,
        "s26_and_action_mismatch_count" => n_overlap_s26_multi,
        "s32_and_action_mismatch_count" => n_overlap_s32_multi,
        "s21_and_s26_violation_overlap_count" => n_overlap_s21_s26,
        "mean_delta_log_b" => mean_logb,
        "sd_delta_log_b" =>
            sqrt(max(0.0,logb_sum_sq/n-mean_logb^2)),
        "min_delta_log_b" => logb_min,
        "max_delta_log_b" => logb_max,
        "min_section21_signed_threshold_change" => threshold_movement_min,
        "max_section21_signed_threshold_change" => threshold_movement_max,
        "three_choice_CES_power_mean" => total_power_3,
        "CES_power_increment_from_mismatched_actions" => extra_power,
        "CES_power_increment_from_section21_violations" => extra_power_on_s21,
        "extra_investment_from_mismatched_actions" => extra_investment,
        "extra_investment_from_section21_violations" => extra_investment_on_s21,
        "extra_investment_from_section26_q_upper_violations" =>
            extra_investment_on_s26,
        "extra_investment_from_section32_violations" =>
            extra_investment_on_s32,
        "delta_phi_from_mismatched_actions" => phi_with_extra-phi3,
        "delta_phi_from_section21_violations_first_order_not_additive" =>
            (total_power_3+extra_power_on_s21)^(inv(rho))-phi3,
        "benchmark_expected_X_minus_three_choice_at_expected_K" =>
            expected_X-at_expected.X,
    )
    denominator = expected_X-at_expected.X
    if abs(denominator)>1e-12
        output["counterfactual_extra_X_divided_by_expected_K_shortfall"] =
            extra_investment/denominator
    end
    return output
end

function main()
    state = initialize_simulation(; N=110000, base_seed=20261101)
    cal = state.calibration
    p = cal.parameters
    s = cal.steady_state
    coefficient = cal.threshold.A^(cal.constants.rho-1.0)/cal.constants.rho
    n = length(state.firms.capital)

    qsum = 0.0
    violated = 0
    for l in state.firms.lambda
        q = -log1p(-p.delta)/log(l)
        qsum += q
        violated += (q > 1.0)
    end
    output = Dict{String,Any}(
        "configuration"=>Dict(
            "firms"=>n, "seed"=>20261101,
            "continuum_K"=>s.K, "continuum_X"=>s.X,
            "multi_step_is_counterfactual" => true,
            "multi_step_model_warning" => "Footnotes 2 and 4 allow finite multi-step adjustment under wider productivity support; NOT the published three-action benchmark.",
            "expected_stationary_investment_participation_ignoring_productivity"=>qsum/n,
            "count_lambda_does_not_satisfy_one_sided_condition"=>violated,
            "equation18_investment_multiplier"=>coefficient,
            "julia_version" => string(VERSION),
            "paper_s26_q_interpretation" =>
                "q from Section 2.4; q<0.5 in that analytical argument; q=0.5 is an upper-envelope diagnostic, not a GE-calibrated q.",
            "normal_support_satisfies_s26_probability_one" => false,
            "lambda_violation_theoretical_probability" =>
                cdf(lambda_distribution(p),inv(1.0-p.delta)),
            "lambda_violation_theoretical_expected_count" =>
                n*cdf(lambda_distribution(p),inv(1.0-p.delta)),
        ),
    )
    periods = Vector{Dict{String,Any}}()
    for t in 0:5
        K_prior = state.K
        expected = compute_expected_equilibrium(
            state.re, cal, K_prior, state.mu_prior
        )
        initial_gap = gap_diagnostics(state)
        draw_productivity!(state.firms, state.rngs, cal)
        enum = enumerate_fixed_points!(state.workspace, state.firms, cal)
        Ke = (1.0-p.delta)*K_prior + coefficient*expected.X_next
        at_expected = response_at_K(state.firms, cal, Ke)
        no_new_productivity = response_at_K(
            state.firms, cal, Ke; use_prior_productivity=true
        )
        support = threshold_support_diagnostics(state, Ke)
        multistep = multi_step_response_at_K(state.firms, cal, Ke)
        assumptions = boundedness_assumption_audit(
            state, Ke, expected.X_next, at_expected, multistep, qsum/n
        )
        info = Dict{String,Any}(
            "t"=>t, "K_prior"=>K_prior,
            "expected_K"=>Ke, "expected_X"=>expected.X_next,
            "reaction_at_expected_K"=>at_expected.Phi,
            "reaction_residual_at_expected_K"=>at_expected.Phi-Ke,
            "investment_at_expected_K"=>at_expected.X,
            "reaction_residual_no_productivity_redraw"=>no_new_productivity.Phi-Ke,
            "investment_no_productivity_redraw"=>no_new_productivity.X,
            "one_step_gap_support_at_expected_K"=>support,
            "paper_assumption_audit_at_expected_K"=>assumptions,
            "multi_step_at_expected_K"=>Dict(
                "reaction_residual"=>multistep.Phi-Ke,
                "aggregate_investment"=>multistep.X,
                "increase_in_reaction_vs_three_action"=>multistep.Phi-at_expected.Phi,
                "increase_in_investment_vs_three_action"=>multistep.X-at_expected.X,
                "n_more_than_one_investment"=>multistep.n_more_than_one_investment,
                "n_more_than_one_divestment"=>multistep.n_more_than_one_divestment,
                "min_exponent"=>multistep.min_exponent,
                "max_exponent"=>multistep.max_exponent,
                "max_gap_violation"=>multistep.max_gap_violation,
            ),
            "investing_at_expected_K"=>at_expected.investing,
            "inactive_at_expected_K"=>at_expected.inactive,
            "divesting_at_expected_K"=>at_expected.divesting,
            "n_fixed"=>enum.total,
            "n_positive_X"=>enum.admissible,
            "min_fixed_K"=>enum.minimum_K, "max_fixed_K"=>enum.maximum_K,
            "min_fixed_X"=>minimum(view(state.workspace.fixed_investment,1:enum.total)),
            "max_fixed_X"=>maximum(view(state.workspace.fixed_investment,1:enum.total)),
            "current_gap"=>initial_gap,
        )
        if enum.admissible == 0
            info["status"] = "no positive investment fixed point; do not change ES1"
            push!(periods, info)
            break
        end
        selection = select_es1(state.workspace,expected.X_next,Ke)
        chosenK = state.workspace.fixed_capital[selection.index]
        chosenX = state.workspace.fixed_investment[selection.index]
        info["status"] = "selected"
        info["K_selected"] = chosenK
        info["X_selected"] = chosenX
        info["epsilon"] = log(chosenX/expected.X_next)
        info["equation18_residual"] = chosenK -
            ((1.0-p.delta)*K_prior + coefficient*chosenX)
        chosen = solve_equilibrium!(
            state.workspace, state.firms, cal, expected.X_next, K_prior
        )
        @assert isapprox(chosen.K_next, chosenK;atol=1e-10,rtol=1e-12)
        state.K = chosen.K_next
        # Existing Issue #5 first-order closure, deliberately unchanged.
        state.mu_prior = s.mu * exp(state.re.h*log(state.K/s.K))
        state.period += 1
        push!(periods, info)
    end
    output["periods"] = periods
    path = normpath(joinpath(@__DIR__,"..","results","initial_collapse_audit.toml"))
    mkpath(dirname(path))
    open(path,"w") do io
        TOML.print(io,output)
    end
    println("Audit saved: ",path)
    println("One-sided-condition violations: ",violated," of ",n)
    for rec in periods
        println("t=",rec["t"]," expected K=",rec["expected_K"],
            " Phi(expected K)-expected K=",rec["reaction_residual_at_expected_K"],
            " selected X=",get(rec,"X_selected","no positive-X fixed point"),
            " fixed points=",rec["n_fixed"])
    end
    return nothing
end

main()
