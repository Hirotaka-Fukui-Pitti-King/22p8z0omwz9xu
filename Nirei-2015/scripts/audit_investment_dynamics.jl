#!/usr/bin/env julia
# Finite-N mechanism audit for the first investment collapse; not a model alternative.
# From Nirei-2015/: julia --project=. scripts/audit_investment_dynamics.jl
using Nirei2015
using TOML
using Distributions: Normal, cdf, ccdf, pdf, truncated
using Random: Xoshiro, rand
using LinearAlgebra: SymTridiagonal, eigen

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


# Four cases at a FIXED candidate aggregate capital.
# A=old thresholds/old CES productivity; B=new thresholds/old CES;
# C=old thresholds/new CES; D=new thresholds/new CES (actual benchmark).
# B and C are identifying counterfactuals, not feasible equilibrium paths.
function four_way_productivity_effects(firms, calibration, K::Float64,
                                       old_response, new_response)
    p = calibration.parameters
    rho = calibration.constants.rho
    n = length(firms.capital)
    power = zeros(Float64,4)
    investment_old = 0.0
    investment_new = 0.0
    investing_old = 0
    investing_new = 0
    @inbounds for i in 1:n
        d = (1.0-p.delta)*firms.capital[i]
        lambda_i = firms.lambda[i]
        old_a = firms.productivity[i]
        new_a = firms.next_productivity[i]
        old_b = threshold_coefficient(lambda_i,old_a,calibration)
        new_b = threshold_coefficient(lambda_i,new_a,calibration)
        old_exponent = K<=d/(lambda_i*old_b) ? -1 :
                       K<=d/old_b ? 0 : 1
        new_exponent = K<=d/(lambda_i*new_b) ? -1 :
                       K<=d/new_b ? 0 : 1
        old_capital = d*lambda_i^old_exponent
        new_capital = d*lambda_i^new_exponent
        old_weight = old_a^(rho/p.alpha)
        new_weight = new_a^(rho/p.alpha)
        power[1] += old_weight*old_capital^rho/n
        power[2] += old_weight*new_capital^rho/n
        power[3] += new_weight*old_capital^rho/n
        power[4] += new_weight*new_capital^rho/n
        investment_old += (old_capital-d)/n
        investment_new += (new_capital-d)/n
        investing_old += old_exponent==1
        investing_new += new_exponent==1
    end
    phi = power.^(inv(rho))
    isapprox(phi[1],old_response.Phi;atol=1e-9,rtol=1e-10) ||
        error("four-way A disagrees with no-redraw reaction")
    isapprox(phi[4],new_response.Phi;atol=1e-9,rtol=1e-10) ||
        error("four-way D disagrees with benchmark reaction")
    isapprox(investment_old,old_response.X;atol=1e-9,rtol=1e-10) ||
        error("four-way no-redraw investment mismatch")
    isapprox(investment_new,new_response.X;atol=1e-9,rtol=1e-10) ||
        error("four-way benchmark investment mismatch")
    threshold_effect = 0.5*((phi[2]-phi[1])+(phi[4]-phi[3]))
    ces_effect = 0.5*((phi[3]-phi[1])+(phi[4]-phi[2]))
    isapprox(threshold_effect+ces_effect,phi[4]-phi[1];
             atol=1e-10,rtol=1e-12) ||
        error("Shapley effects fail to add to total productivity effect")
    return Dict(
        "phi_A_old_threshold_old_CES"=>phi[1],
        "phi_B_new_threshold_old_CES"=>phi[2],
        "phi_C_old_threshold_new_CES"=>phi[3],
        "phi_D_new_threshold_new_CES"=>phi[4],
        "residual_A"=>phi[1]-K,
        "residual_B"=>phi[2]-K,
        "residual_C"=>phi[3]-K,
        "residual_D"=>phi[4]-K,
        "threshold_effect_shapley"=>threshold_effect,
        "CES_weight_effect_shapley"=>ces_effect,
        "total_effect_D_minus_A"=>phi[4]-phi[1],
        "X_A_and_C"=>investment_old,
        "X_B_and_D"=>investment_new,
        "n_investing_old_threshold"=>investing_old,
        "n_investing_new_threshold"=>investing_new,
        "is_fixed_K_counterfactual"=>true,
    )
end

# Analytically integrate over s_0 ~ Uniform[0,1). For each firm,
# k_0(s)=b(a_0,lambda)*Kbar*lambda^s. Next threshold is b(a_1,lambda)*K.
# c=[log(b_1*K/((1-delta)*b_0*Kbar))]/log(lambda).
# Invest s<c; divest s>=c+1; only {-1,0,+1} choices are allowed.
# These are moments conditional on a_0, a_1, lambda and a FIXED K.
# The random finite-N K_0 and hence its conditional expectation are not
# integrated out; this is not an exact ex-ante RE forecast.
@inline function _uniform_gap_integrals(
    lambda_i::Float64, rho::Float64, cutoff::Float64
)
    hi_invest = clamp(cutoff,0.0,1.0)
    lo_divest = clamp(cutoff+1.0,0.0,1.0)
    logl = log(lambda_i)
    g0 = expm1(logl*hi_invest)/logl
    g1 = (lambda_i-exp(logl*lo_divest))/logl
    pow0 = expm1(rho*logl*hi_invest)/(rho*logl)
    pow1 = (exp(rho*logl*lo_divest)-exp(rho*logl*hi_invest))/(rho*logl)
    pow2 = (exp(rho*logl)-exp(rho*logl*lo_divest))/(rho*logl)
    return (
        expected_net_weight=(lambda_i-1.0)*g0+
                            (inv(lambda_i)-1.0)*g1,
        expected_power_weight=lambda_i^rho*pow0+pow1+
                              lambda_i^(-rho)*pow2,
        prob_invest=hi_invest,
        prob_divest=1.0-lo_divest,
    )
end

# Golub-Welsch Gaussian quadrature for Z~N(0,1).
# Recurrence: off-diagonal sqrt(j), j=1,...,n-1.
# No extra dependency and no unverified hard-coded quadrature constants.
function _normal_quadrature(n::Int)
    n>=2 || throw(ArgumentError("normal quadrature requires >=2 nodes"))
    jacobi = SymTridiagonal(zeros(n),sqrt.(Float64.(1:(n-1))))
    eig = eigen(jacobi)
    return eig.values,eig.vectors[1,:].^2
end

# Integrate over initial independent uniform s_0 and next-period iid Normal
# log productivity, conditional on sampled a_0 and lambda. The reference
# prior aggregate capital is Kbar by construction of k_0(s), while
# candidate K is FIXED at the realized conditional forecast.
# In particular, (E[K^rho])^(1/rho) is not E[K].
function initial_uniform_gap_expected_moments(
    firms, calibration, K::Float64; nodes::Int=32
)
    p = calibration.parameters
    rho = calibration.constants.rho
    s = calibration.steady_state
    c = calibration.constants
    gaussian_nodes, gaussian_weights = _normal_quadrature(nodes)
    n = length(firms.capital)
    xsum = 0.0
    powersum = 0.0
    pinvest = 0.0
    pdivest = 0.0
    factor = rho/(p.alpha*(1.0-rho))
    @inbounds for i in 1:n
        lambda_i = firms.lambda[i]
        a0 = firms.productivity[i]
        log_a0 = log(a0)
        b0 = threshold_coefficient(lambda_i,a0,calibration)
        capital_scale = (1.0-p.delta)*b0*s.K
        base_log_ratio = -log1p(-p.delta)+log(K/s.K)
        logl = log(lambda_i)
        for j in eachindex(gaussian_nodes)
            loga1 = c.mean_log_a+p.sigma_log_a*gaussian_nodes[j]
            weight = gaussian_weights[j]/n
            cutoff = (base_log_ratio+factor*(loga1-log_a0))/logl
            moment = _uniform_gap_integrals(lambda_i,rho,cutoff)
            xsum += weight*capital_scale*moment.expected_net_weight
            powersum += weight*exp((rho/p.alpha)*loga1)*
                        capital_scale^rho*moment.expected_power_weight
            pinvest += weight*moment.prob_invest
            pdivest += weight*moment.prob_divest
        end
    end
    return Dict(
        "quadrature_nodes"=>nodes,
        "expected_net_investment"=>xsum,
        "expected_CES_power_mean"=>powersum,
        "CES_capital_of_expected_power_not_expected_capital"=>
            powersum^(inv(rho)),
        "expected_investing_fraction"=>pinvest,
        "expected_divesting_fraction"=>pdivest,
        "conditional_on_initial_lambda_a0_and_fixed_candidate_K"=>true,
        "not_exact_ex_ante_expectation_of_nonlinear_K"=>true,
    )
end


# Continuum reaction for Appendix C's INITIAL uniform gap law (not a
# stationary law of the actual finite-notch Markov process).
#
# With theta=rho/[alpha*(1-rho)] and L_j=log(a_{i,j}), let
# D=L_1-L_0~N(0,2*sigma_a^2).
# E[exp(theta*L_0)|D] =
# exp(theta*(mu-D/2) + theta^2*sigma_a^2/4).
# The CES power weight of next productivity and threshold-based initial
# capital simplifies exactly:
# a_1^(rho/alpha) * a_0^(theta*rho)
# = exp(theta*L_0 + (rho/alpha)*D).
# This collapses the independent old/new productivity double integral
# to a ONE-dimensional Gaussian integral over D.
#
# K0(s)=b(a0,lambda)*Kbar*lambda^s is drawn at the stated initial law.
# The candidate future K is treated as fixed in the reaction. We use
# exactly the same three-action policies as the production equilibrium.
# E[CES power]^(1/rho) is NOT E[finite-N equilibrium K] nor E[ES1 K].
function _legendre_quadrature(n::Int)
    n>=2 || throw(ArgumentError("Legendre quadrature needs >=2 nodes"))
    indices=Float64.(1:(n-1))
    subdiag=indices ./ sqrt.(4.0 .* indices.^2 .- 1.0)
    ev=eigen(SymTridiagonal(zeros(n),subdiag))
    return ev.values, 2.0 .* (ev.vectors[1,:].^2)
end

function _unconditional_initial_continuum_data(
    calibration; lambda_nodes::Int=32, shock_nodes::Int=40
)
    p=calibration.parameters
    distribution=lambda_distribution(p)
    lower=p.lambda_mean-p.lambda_truncation_sd*p.lambda_sd
    upper=p.lambda_mean+p.lambda_truncation_sd*p.lambda_sd
    lambda_z,lambda_w=_legendre_quadrature(lambda_nodes)
    diff_z,diff_w=_normal_quadrature(shock_nodes)
    return (
        lambda_nodes=lower .+ (upper-lower)/2.0 .* (lambda_z .+ 1.0),
        lambda_weights=(upper-lower)/2.0 .* lambda_w .*
            pdf.(Ref(distribution),lower .+
                 (upper-lower)/2.0 .* (lambda_z .+ 1.0)),
        difference_nodes=sqrt(2.0)*p.sigma_log_a .* diff_z,
        difference_weights=diff_w,
        counts=(lambda_nodes,shock_nodes),
    )
end

function unconditional_initial_continuum_reaction(
    calibration, data, K::Float64
)
    p=calibration.parameters
    rho=calibration.constants.rho
    s=calibration.steady_state
    theta=rho/(p.alpha*(1.0-rho))
    K>0.0 || throw(ArgumentError("continuum candidate capital must be positive"))
    ratio=log(K/s.K)
    average_power=0.0
    average_investment=0.0
    investing=0.0
    divesting=0.0
    # Integrate over lambda and difference D of iid Gaussian productivity.
    for k in eachindex(data.lambda_nodes)
        lambda_i=data.lambda_nodes[k]
        weight_lambda=data.lambda_weights[k]
        b0_unit=threshold_coefficient(lambda_i,1.0,calibration)
        scale=(1.0-p.delta)*b0_unit*s.K
        ln_l=log(lambda_i)
        for j in eachindex(data.difference_nodes)
            difference=data.difference_nodes[j]
            conditional_a0_theta=exp(
                theta*(calibration.constants.mean_log_a-0.5*difference)+
                0.25*(theta*p.sigma_log_a)^2
            )
            cutoff=(-log1p(-p.delta)+ratio+theta*difference)/ln_l
            moments=_uniform_gap_integrals(lambda_i,rho,cutoff)
            w=weight_lambda*data.difference_weights[j]*conditional_a0_theta
            average_power+=w*exp((rho/p.alpha)*difference)*
                           scale^rho*moments.expected_power_weight
            average_investment+=w*scale*moments.expected_net_weight
            # Participation fractions must be integrated WITHOUT the
            # conditional capital weight E[a0^theta|D].
            wprob=weight_lambda*data.difference_weights[j]
            investing+=wprob*moments.prob_invest
            divesting+=wprob*moments.prob_divest
        end
    end
    phi=average_power^(inv(rho))
    return (
        Phi=phi, residual=phi-K, X=average_investment,
        investing=investing, divesting=divesting,
        CES_power=average_power,
    )
end

# On a disclosed broad but finite grid, bracket every sign change
# and bisect each bracket. This is evidence for roots in the tested
# interval, NOT a proof that no roots exist elsewhere.
function initial_continuum_consistency_audit(
    calibration, K_expected::Float64;
    lambda_nodes::Int=32, shock_nodes::Int=40,
)
    s=calibration.steady_state
    data=_unconditional_initial_continuum_data(
        calibration;lambda_nodes=lambda_nodes,shock_nodes=shock_nodes
    )
    at_steady=unconditional_initial_continuum_reaction(calibration,data,s.K)
    at_expected=unconditional_initial_continuum_reaction(
        calibration,data,K_expected
    )
    # K is a candidate aggregate within this INITIAL-CROSS-SECTION
    # reaction, not a new steady-state calibration.
    ratios=(0.55,0.65,0.75,0.85,0.90,0.94,0.96,0.975,
            0.985,0.995,1.0,1.005,1.015,1.025,1.04,1.075,1.10,1.20,1.40)
    grid=Vector{Dict{String,Float64}}()
    for ratio in ratios
        v=unconditional_initial_continuum_reaction(
            calibration,data,s.K*ratio
        )
        push!(grid,Dict(
            "K_to_Kbar"=>ratio,"Phi_minus_K"=>v.residual,
            "X"=>v.X,"Phi"=>v.Phi,
        ))
    end
    roots=Float64[]
    root_residuals=Float64[]
    root_investments=Float64[]
    root_investing_fractions=Float64[]
    for z in 1:(length(ratios)-1)
        left=s.K*ratios[z]
        right=s.K*ratios[z+1]
        vl=grid[z]["Phi_minus_K"]
        vr=grid[z+1]["Phi_minus_K"]
        if vl==0.0
            push!(roots,left)
            push!(root_residuals,vl)
            root_stats=unconditional_initial_continuum_reaction(
                calibration,data,left
            )
            push!(root_investments,root_stats.X)
            push!(root_investing_fractions,root_stats.investing)
        elseif (vl<0.0 && vr>0.0) || (vl>0.0 && vr<0.0)
            for iteration in 1:52
                mid=left+(right-left)/2.0
                vm=unconditional_initial_continuum_reaction(
                    calibration,data,mid
                ).residual
                if (vl<0.0 && vm<0.0) || (vl>0.0 && vm>0.0)
                    left=mid
                    vl=vm
                else
                    right=mid
                end
            end
            root=left+(right-left)/2.0
            push!(roots,root)
            root_stats=unconditional_initial_continuum_reaction(
                calibration,data,root
            )
            push!(root_residuals,root_stats.residual)
            push!(root_investments,root_stats.X)
            push!(root_investing_fractions,root_stats.investing)
        end
    end
    return Dict{String,Any}(
        "interpretation"=>"One-period continuum reaction at Appendix C initial distribution, not an invariant stationary law and not an ES1 equilibrium expectation",
        "lambda_quadrature_nodes"=>lambda_nodes,
        "difference_normal_quadrature_nodes"=>shock_nodes,
        "Kbar"=>s.K,
        "K_expected"=>K_expected,
        "Phi_at_Kbar"=>at_steady.Phi,
        "Phi_at_Kbar_minus_Kbar"=>at_steady.residual,
        "X_at_Kbar"=>at_steady.X,
        "X_at_Kbar_minus_Xbar"=>at_steady.X-s.X,
        "investing_at_Kbar"=>at_steady.investing,
        "divesting_at_Kbar"=>at_steady.divesting,
        "Phi_at_expected_K"=>at_expected.Phi,
        "Phi_at_expected_K_minus_expected_K"=>at_expected.residual,
        "X_at_expected_K"=>at_expected.X,
        "X_at_expected_K_minus_expected_Xbar"=>at_expected.X-s.X,
        "root_capitals_in_tested_intervals"=>roots,
        "root_residuals"=>root_residuals,
        "root_investments"=>root_investments,
        "root_investing_fractions"=>root_investing_fractions,
        "tested_ratio_min"=>first(ratios),
        "tested_ratio_max"=>last(ratios),
        "grid"=>grid,
    )
end


# Closed-form Gaussian expectation for the initial one-step continuum.
# Write theta=rho/[alpha*(1-rho)], D=log(a1/a0)~N(0,2*sigma_a^2),
# v(D)=K/[(1-delta)*Kbar]*exp(theta*D), and
# h(D)=clamp(v(D),lambda^-1,lambda).
# Exact s0~Uniform[0,1) integration of ORIGINAL {-1,0,+1} policy:
#  E_s0[X_i | a0,D,lambda]
#     =(1-delta)b0*Kbar * (lambda-1)/log(lambda) * (h(D)-1)
#  E_s0[a1^(rho/alpha)*k1^rho | a0,D,lambda]
#     =a1^(rho/alpha)*[(1-delta)b0*Kbar]^rho *
#       (lambda^rho-1)/(rho*log(lambda)) * h(D)^rho.
# Since E[a0^theta | D] = exp(theta*(mu-D/2) + theta^2*sigma_a^2/4),
# integrate D using tilted Gaussian truncated moments, with NO
# approximation to the Gaussian tails or to clipping.
function _tilted_normal_interval(
    sd::Float64, tilt::Float64, low::Float64, high::Float64
)
    normal=Normal()
    shift=tilt*sd^2
    p_low=isfinite(low) ? cdf(normal,(low-shift)/sd) : 0.0
    p_high=isfinite(high) ? cdf(normal,(high-shift)/sd) : 1.0
    return exp(0.5*(tilt*sd)^2)*max(0.0,p_high-p_low)
end

function _bounded_exponential_normal_moment(
    sd::Float64, tilt::Float64, exponent::Float64,
    base::Float64, theta::Float64, ln_lambda::Float64,
)
    threshold_low=(-ln_lambda-base)/theta
    threshold_high=(ln_lambda-base)/theta
    return exp(-exponent*ln_lambda)*
             _tilted_normal_interval(sd,tilt,-Inf,threshold_low) +
           exp(exponent*base)*
             _tilted_normal_interval(
                 sd,tilt+exponent*theta,threshold_low,threshold_high
             ) +
           exp(exponent*ln_lambda)*
             _tilted_normal_interval(sd,tilt,threshold_high,Inf)
end

function exact_gaussian_initial_continuum_reaction(
    calibration, data, K::Float64
)
    p=calibration.parameters
    rho=calibration.constants.rho
    s=calibration.steady_state
    K>0.0 || throw(ArgumentError("candidate K must be positive"))
    theta=rho/(p.alpha*(1.0-rho))
    sd=sqrt(2.0)*p.sigma_log_a
    base=-log1p(-p.delta)+log(K/s.K)
    prefactor=exp(theta*calibration.constants.mean_log_a+
                  0.25*(theta*p.sigma_log_a)^2)
    power_total=0.0
    x_total=0.0
    for j in eachindex(data.lambda_nodes)
        lambda_i=data.lambda_nodes[j]
        w=data.lambda_weights[j]
        logl=log(lambda_i)
        b0_unit=threshold_coefficient(lambda_i,1.0,calibration)
        scale=(1.0-p.delta)*b0_unit*s.K
        investment_tilt=-0.5*theta
        power_tilt=(rho/p.alpha)-0.5*theta
        exph=_bounded_exponential_normal_moment(
            sd,investment_tilt,1.0,base,theta,logl
        )
        baseline=exp(0.5*(investment_tilt*sd)^2)
        exph_rho=_bounded_exponential_normal_moment(
            sd,power_tilt,rho,base,theta,logl
        )
        x_total+=w*prefactor*scale*(lambda_i-1.0)/logl*
                 (exph-baseline)
        power_total+=w*prefactor*scale^rho*
                     expm1(rho*logl)/(rho*logl)*exph_rho
    end
    phi=power_total^(inv(rho))
    return (Phi=phi,residual=phi-K,X=x_total,CES_power=power_total)
end

# Bracket and bisect the analytic-in-productivity continuum map over the
# same explicitly disclosed range as the independent 2D quadrature.
function exact_initial_continuum_consistency_audit(
    calibration, K_expected::Float64; lambda_nodes::Int=64
)
    s=calibration.steady_state
    data=_unconditional_initial_continuum_data(
        calibration;lambda_nodes=lambda_nodes,shock_nodes=2
    )
    at_bar=exact_gaussian_initial_continuum_reaction(calibration,data,s.K)
    at_exp=exact_gaussian_initial_continuum_reaction(
        calibration,data,K_expected
    )
    ratios=(0.55,0.65,0.75,0.85,0.90,0.94,0.96,0.975,
            0.985,0.995,1.0,1.005,1.015,1.025,1.04,1.075,1.10,1.20,1.40)
    residuals=Float64[]
    roots=Float64[]
    x_roots=Float64[]
    for ratio in ratios
        push!(residuals,
            exact_gaussian_initial_continuum_reaction(
                calibration,data,s.K*ratio
            ).residual)
    end
    for j in 1:(length(ratios)-1)
        lo=s.K*ratios[j]
        hi=s.K*ratios[j+1]
        flo=residuals[j]
        fhi=residuals[j+1]
        if flo==0.0
            push!(roots,lo)
            push!(x_roots,
                exact_gaussian_initial_continuum_reaction(
                    calibration,data,lo
                ).X)
        elseif (flo<0.0 && fhi>0.0)||(flo>0.0 && fhi<0.0)
            for iteration in 1:52
                middle=lo+0.5*(hi-lo)
                fmid=exact_gaussian_initial_continuum_reaction(
                    calibration,data,middle
                ).residual
                if (flo<0.0 && fmid<0.0)||(flo>0.0 && fmid>0.0)
                    lo=middle
                    flo=fmid
                else
                    hi=middle
                end
            end
            root=lo+0.5*(hi-lo)
            push!(roots,root)
            push!(x_roots,
                exact_gaussian_initial_continuum_reaction(
                    calibration,data,root
                ).X)
        end
    end
    return Dict{String,Any}(
        "method"=>"exact Gaussian difference integration with normal CDF, Gauss-Legendre only for lambda",
        "Kbar"=>s.K,
        "lambda_nodes"=>lambda_nodes,
        "Phi_at_Kbar_minus_Kbar"=>at_bar.residual,
        "X_at_Kbar_minus_Xbar"=>at_bar.X-s.X,
        "Phi_at_expected_K_minus_expected_K"=>at_exp.residual,
        "X_at_expected_K"=>at_exp.X,
        "root_capitals_in_tested_intervals"=>roots,
        "root_investments"=>x_roots,
        "searched_K_to_Kbar_range_min"=>first(ratios),
        "searched_K_to_Kbar_range_max"=>last(ratios),
        "not_a_long_run_stationary_distribution"=>true,
    )
end


# Distribution-transition diagnostic. The candidate aggregate capital is
# held exogenous at published Kbar; neither a finite-N GE equilibrium nor
# an alternative ES1 trajectory is constructed. All firm capital updates
# follow the ORIGINAL {-1,0,+1} choice and iid productivity redraw law.
#
# This tests whether the INITIAL joint law of actual (unwrapped) gaps
# survives transitions, and whether observed fixed-K block moments are
# approximately stable. It does NOT prove existence/uniqueness of a
# self-consistent invariant distribution for endogenous aggregate capital.
function _fixed_K_distribution_snapshot(
    firms, calibration, Kfixed::Float64, date::Int,
    net_investment::Float64, capital_identity_error::Float64,
    investing_fraction::Float64, divesting_fraction::Float64,
)
    n=length(firms.capital)
    rho=calibration.constants.rho
    alpha=calibration.parameters.alpha
    sum_k=0.0
    sum_power=0.0
    sum_s=0.0
    sum_s_squared=0.0
    outside=0
    below=0
    above=0
    bins=zeros(Int,10)
    sum_loga=0.0
    sum_logk=0.0
    sum_loga2=0.0
    sum_logk2=0.0
    sum_logak=0.0
    @inbounds for i in 1:n
        k_i=firms.capital[i]
        a_i=firms.productivity[i]
        lambda_i=firms.lambda[i]
        b_i=threshold_coefficient(lambda_i,a_i,calibration)
        actual_gap=log(k_i/(b_i*Kfixed))/log(lambda_i)
        fractional_gap=mod(actual_gap,1.0)
        bin_idx=clamp(floor(Int,fractional_gap*10.0)+1,1,10)
        bins[bin_idx]+=1
        outside+=(actual_gap<0.0 || actual_gap>=1.0)
        below+=(actual_gap<0.0)
        above+=(actual_gap>=1.0)
        sum_s+=actual_gap
        sum_s_squared+=actual_gap^2
        sum_k+=k_i
        sum_power+=a_i^(rho/alpha)*k_i^rho
        la=log(a_i)
        lk=log(k_i)
        sum_loga+=la
        sum_logk+=lk
        sum_loga2+=la^2
        sum_logk2+=lk^2
        sum_logak+=la*lk
    end
    mean_s=sum_s/n
    mean_la=sum_loga/n
    mean_lk=sum_logk/n
    var_la=max(0.0,sum_loga2/n-mean_la^2)
    var_lk=max(0.0,sum_logk2/n-mean_lk^2)
    correlation=(var_la>0.0 && var_lk>0.0) ?
        (sum_logak/n-mean_la*mean_lk)/sqrt(var_la*var_lk) : 0.0
    return Dict{String,Any}(
        "t"=>date,
        "mean_firm_capital"=>sum_k/n,
        "CES_capital"=>((sum_power/n)^(inv(rho))),
        "CES_minus_exogenous_Kbar"=>(sum_power/n)^(inv(rho))-Kfixed,
        "gross_investment_net_of_depreciated_capital"=>net_investment,
        "capital_identity_error"=>capital_identity_error,
        "fraction_investing"=>investing_fraction,
        "fraction_divesting"=>divesting_fraction,
        "fraction_actual_gap_outside_0_1"=>outside/n,
        "fraction_actual_gap_below_zero"=>below/n,
        "fraction_actual_gap_at_least_one"=>above/n,
        "mean_actual_gap"=>mean_s,
        "sd_actual_gap"=>sqrt(max(0.0,sum_s_squared/n-mean_s^2)),
        "fractional_gap_bin_frequencies"=>bins ./ n,
        "max_fractional_gap_bin_error_from_uniform"=>
            maximum(abs.(bins ./ n .- 0.1)),
        "correlation_log_productivity_log_firm_capital"=>correlation,
    )
end


# One EXACT firm transition under an exogenously fixed aggregate K. Shared
# by the earlier distribution audit and the new initial-distribution
# preconditioning experiment. Not a GE solver or an equilibrium selection.
function _fixed_K_firm_transition!(
    firms, rngs, calibration, Kfixed::Float64
)
    n=length(firms.capital)
    p=calibration.parameters
    draw_productivity!(firms,rngs,calibration)
    sum_old=0.0
    sum_new=0.0
    total_x=0.0
    investing=0
    divesting=0
    @inbounds for i in eachindex(firms.capital)
        k_old=firms.capital[i]
        d=(1.0-p.delta)*k_old
        lambda_i=firms.lambda[i]
        a_new=firms.next_productivity[i]
        b_next=threshold_coefficient(lambda_i,a_new,calibration)
        exponent=Kfixed<=d/(lambda_i*b_next) ? -1 :
                 Kfixed<=d/b_next ? 0 : 1
        k_new=d*lambda_i^exponent
        firms.capital[i]=k_new
        sum_old+=k_old
        sum_new+=k_new
        total_x+=k_new-d
        investing+=(exponent==1)
        divesting+=(exponent==-1)
    end
    advance_productivity!(firms)
    old_mean=sum_old/n
    new_mean=sum_new/n
    xmean=total_x/n
    identity_error=new_mean-(1.0-p.delta)*old_mean-xmean
    return (
        old_mean=old_mean,new_mean=new_mean,
        xmean=xmean,identity_error=identity_error,
        investing=investing,divesting=divesting,
    )
end

function fixed_K_distribution_diagnostic(
    calibration; n::Int=110000, seed::Int=20261101, horizon::Int=256
)
    horizon>=128 || throw(ArgumentError("fixed-K horizon must be >=128"))
    rngs=firm_rngs(seed,1)
    firms=initialize_firms(calibration,rngs;N=n)
    p=calibration.parameters
    Kfixed=calibration.steady_state.K
    checks=(0,1,2,5,10,25,50,100,128,192,256)
    snapshots=Vector{Dict{String,Any}}()
    push!(snapshots,_fixed_K_distribution_snapshot(
        firms,calibration,Kfixed,0,0.0,0.0,0.0,0.0
    ))
    block_x=zeros(Float64,2)
    block_capital=zeros(Float64,2)
    block_net_growth=zeros(Float64,2)
    block_count=zeros(Int,2)
    largest_identity_error=0.0
    for t in 1:horizon
        step=_fixed_K_firm_transition!(firms,rngs,calibration,Kfixed)
        old_mean=step.old_mean
        new_mean=step.new_mean
        xmean=step.xmean
        identity_error=step.identity_error
        investing=step.investing
        divesting=step.divesting
        largest_identity_error=max(largest_identity_error,abs(identity_error))
        if t>horizon-128
            block=t<=horizon-64 ? 1 : 2
            block_x[block]+=xmean
            block_capital[block]+=new_mean
            block_net_growth[block]+=new_mean-old_mean
            block_count[block]+=1
        end
        if t in checks
            push!(snapshots,_fixed_K_distribution_snapshot(
                firms,calibration,Kfixed,t,xmean,identity_error,
                investing/n,divesting/n
            ))
        end
    end
    block_x ./= block_count
    block_capital ./= block_count
    block_net_growth ./= block_count
    return Dict{String,Any}(
        "interpretation"=>"Diagnostic exogenous K=Kbar, original 3-choice and iid shocks; not equilibrium, not ES1, not proof of stationarity",
        "seed"=>seed,
        "firms"=>n,
        "horizon"=>horizon,
        "K_fixed_exogenous"=>Kfixed,
        "snapshots"=>snapshots,
        "last_two_64_period_mean_gross_investment"=>block_x,
        "last_two_64_period_mean_firm_capital"=>block_capital,
        "last_two_64_period_mean_change_in_firm_capital"=>block_net_growth,
        "last_block_mean_gross_X_minus_delta_mean_k_approx"=>
            block_x[2]-p.delta*block_capital[2],
        "largest_capital_accounting_identity_error"=>largest_identity_error,
    )
end


# Case B is a fully labeled SUPPORT-ONLY counterfactual: the published
# continuum calibration (A, Kbar, linear RE coefficients) remains identical
# to Case A. Re-solving a new continuum calibration would be a SECOND change.
# Conditional replacement coupling has the correct independent conditional
# lambda marginal: accept each original lambda when above cutoff; otherwise
# independently redraw from the original law conditional on lambda>cutoff.
# All a0/s0 and subsequent productivity RNG streams are exactly unchanged.
function _apply_lambda_support_counterfactual!(sim_state)
    cal=sim_state.calibration
    p=cal.parameters
    s=cal.steady_state
    floor_value=nextfloat(inv(1.0-p.delta))
    upper=p.lambda_mean+p.lambda_truncation_sd*p.lambda_sd
    original_lower=p.lambda_mean-p.lambda_truncation_sd*p.lambda_sd
    original_lower<floor_value<upper ||
        throw(ArgumentError("support-only test requires the cutoff to be inside the published support"))
    corrected_distribution=truncated(
        Normal(p.lambda_mean,p.lambda_sd),floor_value,upper
    )
    replacement_rng=Xoshiro(UInt64(sim_state.base_seed) ⊻ 0x397e1b3a935e57af)
    changed=Int[]
    for i in eachindex(sim_state.firms.lambda)
        if sim_state.firms.lambda[i]<=inv(1.0-p.delta)
            lambda_new=rand(replacement_rng,corrected_distribution)
            lambda_new>inv(1.0-p.delta) ||
                error("The conditional lambda sampler hit the excluded boundary")
            sim_state.firms.lambda[i]=lambda_new
            b=threshold_coefficient(
                lambda_new,sim_state.firms.productivity[i],cal
            )
            sim_state.firms.initial_b[i]=b
            sim_state.firms.capital[i]=
                b*s.K*exp(sim_state.firms.initial_gap[i]*log(lambda_new))
            push!(changed,i)
        end
    end
    all(l->l*(1.0-p.delta)>1.0,sim_state.firms.lambda) ||
        error("Case B still violates lambda*(1-delta)>1")
    sim_state.K=aggregate_capital(sim_state.firms,cal)
    sim_state.mu_prior=s.mu*exp(sim_state.re.h*log(sim_state.K/s.K))
    return changed
end

function _support_sample_description(state,modified_count::Int)
    l=state.firms.lambda
    n=length(l)
    avg=sum(l)/n
    sd=sqrt(sum((v-avg)^2 for v in l)/(n-1))
    return Dict{String,Any}(
        "mean"=>avg,"sd"=>sd,"minimum"=>minimum(l),
        "maximum"=>maximum(l),"modified_firms"=>modified_count,
        "violations"=>count(v->v*(1.0-state.calibration.parameters.delta)<=1.0,l)
    )
end

# The comparison is at the SAME candidate K, for the EXACT subset of
# firms replaced in Case B. The other firms have identical capital,
# lambda, shocks and policies. This is the direct (fixed-K) channel, NOT
# a decomposition of the generally discontinuous ES1-selected outcome.
function _direct_changed_firm_effect(stateA,stateB,indices,Ka::Float64)
    cal=stateA.calibration
    p=cal.parameters
    rho=cal.constants.rho
    n=length(stateA.firms.lambda)
    change_x=0.0
    change_power=0.0
    modified_indices=Set(indices)
    untouched_count=0
    for i in 1:n
        if !(i in modified_indices)
            stateA.firms.lambda[i]==stateB.firms.lambda[i] ||
                error("Unchanged lambda accidentally modified")
            stateA.firms.capital[i]==stateB.firms.capital[i] ||
                error("Unchanged firm capital accidentally modified")
            untouched_count+=1
            continue
        end
        a=stateA.firms.next_productivity[i]
        a==stateB.firms.next_productivity[i] ||
            error("Common productivity shocks not preserved")
        x=Vector{Float64}(undef,2)
        power=Vector{Float64}(undef,2)
        for (j,st) in enumerate((stateA,stateB))
            lambda_i=st.firms.lambda[i]
            d=(1.0-p.delta)*st.firms.capital[i]
            b=threshold_coefficient(lambda_i,a,cal)
            adjustment=Ka<=d/(lambda_i*b) ? -1 :
                       Ka<=d/b ? 0 : 1
            knext=d*lambda_i^adjustment
            x[j]=knext-d
            power[j]=a^(rho/p.alpha)*knext^rho
        end
        change_x+=(x[2]-x[1])/n
        change_power+=(power[2]-power[1])/n
    end
    a_response=response_at_K(stateA.firms,cal,Ka)
    b_response=response_at_K(stateB.firms,cal,Ka)
    isapprox(b_response.X-a_response.X,change_x;atol=1e-10,rtol=1e-8) ||
        error("Changed-firm investment attribution does not reconcile")
    isapprox(b_response.Phi^rho-a_response.Phi^rho,change_power;
             atol=1e-10,rtol=1e-8) ||
        error("Changed-firm CES-power attribution does not reconcile")
    return Dict(
        "fixed_candidate_K"=>Ka,
        "changed_firms"=>length(indices),
        "unchanged_firms"=>untouched_count,
        "X_direct_change_at_fixed_K"=>change_x,
        "CES_power_direct_change_at_fixed_K"=>change_power,
        "Phi_direct_change_at_fixed_K"=>b_response.Phi-a_response.Phi,
        "all_unchanged_firm_contributions_identical"=>true,
        "not_selected_equilibrium_attribution"=>true,
    )
end

function _first_equilibrium_support_summary(state)
    cal=state.calibration
    expected=compute_expected_equilibrium(
        state.re,cal,state.K,state.mu_prior
    )
    p=cal.parameters
    Ke=(1.0-p.delta)*state.K+
       cal.threshold.A^(cal.constants.rho-1.0)/
       cal.constants.rho*expected.X_next
    draw_productivity!(state.firms,state.rngs,cal)
    enumeration=enumerate_fixed_points!(
        state.workspace,state.firms,cal
    )
    at_expected=response_at_K(state.firms,cal,Ke)
    result=Dict{String,Any}(
        "K0"=>state.K,"expected_X"=>expected.X_next,
        "expected_K"=>Ke,
        "reaction_residual_at_expected_K"=>at_expected.Phi-Ke,
        "investment_at_expected_K"=>at_expected.X,
        "fixed_points"=>enumeration.total,
        "positive_investment_fixed_points"=>enumeration.admissible,
        "n_investing_at_expected_K"=>at_expected.investing,
        "n_divesting_at_expected_K"=>at_expected.divesting,
    )
    if enumeration.admissible>0
        selected=select_es1(
            state.workspace,expected.X_next,Ke
        )
        result["ES1_X"]=state.workspace.fixed_investment[selected.index]
        result["ES1_K"]=state.workspace.fixed_capital[selected.index]
    else
        result["ES1_status"]="No positive-X fixed point"
    end
    return result
end

function _support_case_path(;horizon::Int=20)
    sim_state=initialize_simulation(;N=110000,base_seed=20261101)
    changed=_apply_lambda_support_counterfactual!(sim_state)
    series=Vector{Dict{String,Any}}()
    failure="none"
    for t in 0:(horizon-1)
        try
            rec=simulate_period!(sim_state)
            push!(series,Dict{String,Any}(
                "t"=>t,"K"=>rec.K_next_firm,
                "X"=>rec.X,"expected_X"=>rec.expected_X,
                "expected_K"=>rec.expected_K_next,
                "epsilon"=>rec.epsilon,"n_fixed"=>rec.n_fixed,
                "n_positive_X"=>rec.n_positive_X,
                "n_investing"=>rec.n_investing,
                "n_divesting"=>rec.n_divesting,
                "fixed_point_residual"=>rec.fixed_point_residual,
                "investment_residual"=>rec.investment_residual,
            ))
        catch error_object
            failure=sprint(showerror,error_object)
            break
        end
    end
    return Dict{String,Any}(
        "requested_periods"=>horizon,
        "completed_periods"=>length(series),
        "completed"=>length(series)==horizon,
        "failure"=>failure,
        "periods"=>series,
        "changed_firms"=>length(changed),
    )
end

function support_consistent_comparison(;horizon::Int=20)
    a=initialize_simulation(;N=110000,base_seed=20261101)
    b=initialize_simulation(;N=110000,base_seed=20261101)
    a.firms.productivity==b.firms.productivity ||
        error("Initial productivity random numbers differ")
    a.firms.initial_gap==b.firms.initial_gap ||
        error("Initial gap random numbers differ")
    sampleA=_support_sample_description(a,0)
    changed=_apply_lambda_support_counterfactual!(b)
    sampleB=_support_sample_description(b,length(changed))
    sampleA["violations"]==673 ||
        error("Unexpected benchmark support violation count: "*
              string(sampleA["violations"]))
    sampleB["violations"]==0 || error("Case B violates support")
    firstA=_first_equilibrium_support_summary(a)
    firstB=_first_equilibrium_support_summary(b)
    a.firms.next_productivity==b.firms.next_productivity ||
        error("Subsequent productivity draws differ in first period")
    isapprox(firstA["ES1_X"],0.06936827009543808;
             atol=1e-9,rtol=1e-10) ||
        error("Baseline Case A first-period ES1 investment regression")
    isapprox(firstA["ES1_K"],8.546731346350889;
             atol=1e-9,rtol=1e-10) ||
        error("Baseline Case A first-period ES1 capital regression")
    direct=_direct_changed_firm_effect(
        a,b,changed,firstA["expected_K"]
    )
    path=_support_case_path(;horizon=horizon)
    return Dict{String,Any}(
        "interpretation"=>"Case A published calibration vs conditional replacement of violating lambdas only; macro calibration held fixed in both; diagnostic Case B NOT published Appendix C",
        "case_A_distribution"=>sampleA,
        "case_B_distribution"=>sampleB,
        "case_A_initial_equilibrium"=>firstA,
        "case_B_initial_equilibrium"=>firstB,
        "direct_modified_firm_effect_at_A_expected_K"=>direct,
        "case_B_simulation"=>path,
        "initial_productivity_shocks_identical"=>true,
        "initial_gaps_identical"=>true,
        "all_subsequent_productivity_streams_identical_by_construction"=>true,
        "published_calibration_unchanged"=>true,
        "case_B_uses_theory_support_but_still_unbounded_productivity"=>true,
    )
end


# PRECONDITIONED-DISTRIBUTION COUNTERFACTUAL. The initial firms are
# exactly the benchmark initial firms; we apply the already audited
# fixed-K firm transition for 256 dates using a SEPARATE shock RNG.
# Original future GE shock RNG state remains at its initial draw.
# Case A and B therefore receive identical a_{i,t+1} for every common
# post-initialization simulation date. Neither this function nor the
# baseline implementation repositions gaps or rescales capital.
function _precondition_simulation_state!(sim_state;warmup::Int=256)
    warmup>=1 || throw(ArgumentError("positive warmup required"))
    base_seed=sim_state.base_seed
    pre_rngs=firm_rngs(base_seed+100000003,sim_state.replication)
    cal=sim_state.calibration
    Kfixed=cal.steady_state.K
    capital_before=aggregate_capital(sim_state.firms,cal)
    initial_lambda=copy(sim_state.firms.lambda)
    last_step=nothing
    for _ in 1:warmup
        last_step=_fixed_K_firm_transition!(
            sim_state.firms,pre_rngs,cal,Kfixed
        )
    end
    sim_state.firms.lambda==initial_lambda ||
        error("Firm indivisibilities changed during preconditioning")
    sim_state.K=aggregate_capital(sim_state.firms,cal)
    sim_state.mu_prior=cal.steady_state.mu*
        exp(sim_state.re.h*log(sim_state.K/cal.steady_state.K))
    sim_state.period=0
    return Dict{String,Any}(
        "warmup_periods"=>warmup,
        "separate_preconditioning_seed"=>base_seed+100000003,
        "exogenous_K_during_warmup"=>Kfixed,
        "CES_capital_before"=>capital_before,
        "CES_capital_after"=>sim_state.K,
        "firm_capital_mean_after"=>
            sum(sim_state.firms.capital)/length(sim_state.firms.capital),
        "last_preconditioning_step_X"=>last_step.xmean,
        "last_preconditioning_step_capital_identity_error"=>
            last_step.identity_error,
        "future_GE_rng_unchanged"=>true,
        "initial_mu_on_original_linear_stable_eigenspace"=>true,
    )
end

function _preconditioned_first_period_diagnostic()
    caseA=initialize_simulation(;N=110000,base_seed=20261101)
    caseB=initialize_simulation(;N=110000,base_seed=20261101)
    caseA.firms.capital==caseB.firms.capital ||
        error("Case A and B do not start with identical firm capital")
    caseA.firms.productivity==caseB.firms.productivity ||
        error("Case A and B initial productivity differ")
    caseA.firms.initial_gap==caseB.firms.initial_gap ||
        error("Case A and B initial gap differ")
    warmup=_precondition_simulation_state!(caseB;warmup=256)
    beforeA=gap_diagnostics(caseA)
    beforeB=gap_diagnostics(caseB)
    firstA=_first_equilibrium_support_summary(caseA)
    firstB=_first_equilibrium_support_summary(caseB)
    caseA.firms.next_productivity==caseB.firms.next_productivity ||
        error("Common post-warmup t=0 productivity innovations DIFFER")
    firstA["ES1_X"]==0.06936827009543808 ||
        error("Published Case A first-period investment changed")
    firstA["ES1_K"]==8.546731346350889 ||
        error("Published Case A first-period capital changed")
    return Dict{String,Any}(
        "caseA"=>firstA,
        "caseB"=>firstB,
        "caseA_mean_firm_capital"=>
            sum(caseA.firms.capital)/length(caseA.firms.capital),
        "caseB_mean_firm_capital"=>
            sum(caseB.firms.capital)/length(caseB.firms.capital),
        "caseA_initial_gap"=>beforeA,
        "caseB_preconditioned_gap"=>beforeB,
        "preconditioning"=>warmup,
        "common_post_GE_shock_vectors_at_t0"=>true,
        "caseB_is_not_published_initialization"=>true,
    )
end

# In contrast to the one-period summary above, this invokes the ORIGINAL
# simulate_period! sequence without disabling macro closure, ES1, or
# positive-investment requirement. Report the real stopping date.
function _preconditioned_GE_path(;horizon::Int=20)
    sim_state=initialize_simulation(;N=110000,base_seed=20261101)
    warmup=_precondition_simulation_state!(sim_state;warmup=256)
    records=Vector{Dict{String,Any}}()
    failure="none"
    for _ in 1:horizon
        try
            rec=simulate_period!(sim_state)
            push!(records,Dict{String,Any}(
                "t"=>rec.t,
                "K"=>rec.K_next_firm,
                "X"=>rec.X,
                "expected_X"=>rec.expected_X,
                "expected_K"=>rec.expected_K_next,
                "epsilon"=>rec.epsilon,
                "n_fixed"=>rec.n_fixed,
                "n_positive_X"=>rec.n_positive_X,
                "n_investing"=>rec.n_investing,
                "n_divesting"=>rec.n_divesting,
                "fixed_point_residual"=>rec.fixed_point_residual,
                "investment_residual"=>rec.investment_residual,
                "logK_expectation_error"=>rec.logK_expectation_error,
            ))
        catch caught
            failure=sprint(showerror,caught)
            break
        end
    end
    return Dict{String,Any}(
        "interpretation"=>"Diagnostic preconditioning at exogenous Kbar; not the published simulation initialization",
        "preconditioning"=>warmup,
        "requested_periods"=>horizon,
        "completed_periods"=>length(records),
        "completed"=>length(records)==horizon,
        "failure"=>failure,
        "periods"=>records,
        "final_GE_capital"=>sim_state.K,
        "not_published_replication"=>true,
    )
end

function preconditioned_distribution_comparison()
    first_period=_preconditioned_first_period_diagnostic()
    short_run=_preconditioned_GE_path(;horizon=20)
    result=Dict{String,Any}(
        "first_period"=>first_period,
        "caseB_twenty_periods"=>short_run,
        "caseA_known_baseline_failure_t"=>5,
        "common_future_productivity_streams_by_separate_precondition_rng"=>true,
    )
    if short_run["completed"]
        result["caseB_1100_periods"]=_preconditioned_GE_path(
            ;horizon=1100
        )
    else
        result["caseB_1100_status"]=
            "Not attempted because preconditioned 20-period gate failed"
    end
    return result
end


# Issue #5: EQUILIBRIUM-SELECTION PATH DEPENDENCE.
# All alternative rules are strictly DIAGNOSTIC. The production ES1
# in src/equilibrium.jl and simulation state recursion are UNCHANGED.
# The fixed-point list is ordered by capital (a strictly monotone
# sweep). The investment list is weakly monotone for the same reason.
# Save every candidate of the original ES1 benchmark for auditability.
function _selection_ranking(workspace, expected_X::Float64, expected_K::Float64)
    nfixed=workspace.n_fixed
    nfixed>0 || error("Fixed points must be enumerated")
    capital=workspace.fixed_capital
    investment=workspace.fixed_investment
    # Exact log distance, all candidates, as in Definition 1.
    best_cap=1
    best_cap_score=abs(log(capital[1]/expected_K))
    for j in 2:nfixed
        sc=abs(log(capital[j]/expected_K))
        if sc<best_cap_score ||
           (sc==best_cap_score && capital[j]<capital[best_cap])
            best_cap=j
            best_cap_score=sc
        end
    end
    # Section 3.4 selection and deterministic second-best positive X.
    valid=Int[]
    for j in 1:nfixed
        if investment[j]>0.0
            push!(valid,j)
        end
    end
    sort!(valid;by=j->(
        abs(log(investment[j]/expected_X)),
        abs(log(capital[j]/expected_K)),capital[j]
    ))
    invest_best=isempty(valid) ? 0 : valid[1]
    second_best=length(valid)>=2 ? valid[2] : 0
    max_positive=isempty(valid) ? 0 : maximum(valid)
    return (
        capital_best=best_cap,
        invest_best=invest_best,
        max_positive=max_positive,
        second_best=second_best,
        positive=length(valid),
        nfixed=nfixed
    )
end

function _all_fixed_points_ranked(workspace, expected_X::Float64, expected_K::Float64)
    r=_selection_ranking(workspace,expected_X,expected_K)
    entries=Vector{Dict{String,Any}}(undef,r.nfixed)
    for j in 1:r.nfixed
        k=workspace.fixed_capital[j]
        x=workspace.fixed_investment[j]
        entries[j]=Dict{String,Any}(
            "capital_rank_ascending"=>j,
            "capital"=>k,"investment"=>x,
            "positive_investment"=>x>0.0,
            "absolute_logK_gap"=>abs(log(k/expected_K)),
            "absolute_logX_gap_or_NaN"=>x>0.0 ?
                abs(log(x/expected_X)) : "undefined for nonpositive X",
        )
    end
    return entries
end

function _period_selection_audit(state)
    cal=state.calibration
    p=cal.parameters
    oldK=state.K
    expected=compute_expected_equilibrium(
        state.re,cal,oldK,state.mu_prior
    )
    expected_X=expected.X_next
    expected_K=(1.0-p.delta)*oldK+
        cal.threshold.A^(cal.constants.rho-1.0)/
        cal.constants.rho*expected_X
    draw_productivity!(state.firms,state.rngs,cal)
    enumeration=enumerate_fixed_points!(state.workspace,state.firms,cal)
    ranks=_selection_ranking(state.workspace,expected_X,expected_K)
    ws=state.workspace
    ranks.positive==enumeration.admissible ||
        error("Positive fixed point count mismatch")
    if ranks.invest_best!=0
        current=select_es1(ws,expected_X,expected_K)
        current.index==ranks.invest_best ||
            error("Diagnostic ES1 disagrees with published selector")
    end
    record=Dict{String,Any}(
        "date"=>state.period,
        "K_prior"=>oldK,
        "expected_X"=>expected_X,
        "expected_K"=>expected_K,
        "total_fixed_points"=>ranks.nfixed,
        "positive_investment_fixed_points"=>ranks.positive,
        "minimum_K"=>minimum(view(ws.fixed_capital,1:ranks.nfixed)),
        "maximum_K"=>maximum(view(ws.fixed_capital,1:ranks.nfixed)),
        "minimum_X"=>minimum(view(ws.fixed_investment,1:ranks.nfixed)),
        "maximum_X"=>maximum(view(ws.fixed_investment,1:ranks.nfixed)),
        "capital_distance_index"=>ranks.capital_best,
        "capital_distance_X"=>ws.fixed_investment[ranks.capital_best],
        "capital_distance_K"=>ws.fixed_capital[ranks.capital_best],
        "capital_and_investment_rules_agree"=>
            ranks.invest_best==ranks.capital_best,
        "all_fixed_points"=>_all_fixed_points_ranked(ws,expected_X,expected_K),
    )
    if ranks.invest_best>0
        record["investment_ES1_index"]=ranks.invest_best
        record["investment_ES1_capital_rank_ascending"]=ranks.invest_best
        record["investment_ES1_investment_rank_ascending"]=ranks.invest_best
        record["investment_ES1_K"]=ws.fixed_capital[ranks.invest_best]
        record["investment_ES1_X"]=ws.fixed_investment[ranks.invest_best]
        record["max_capital_positive_index"]=ranks.max_positive
        record["second_best_positive_index"]=ranks.second_best
        record["ES1_is_max_capital"]=ranks.invest_best==ranks.nfixed
        record["ES1_is_max_positive_capital"]=
            ranks.invest_best==ranks.max_positive
    else
        record["investment_ES1_status"]="no admissible X>0"
    end
    return record,ranks,expected_X,expected_K
end

# Reuses the EXACT candidate fixed points from production enumeration.
# The selected candidate can be any index for DIAGNOSTIC purposes,
# but must have X>0 to retain the Section 3.4 log-shock definition.
# Capital is recomputed independently with eq(13); X with eq(14).
function _advance_with_diagnostic_selection!(
    state, chosen_index::Int, expected_X::Float64
)
    ws=state.workspace
    cal=state.calibration
    p=cal.parameters
    rho=cal.constants.rho
    n=length(state.firms.capital)
    1<=chosen_index<=ws.n_fixed ||
        throw(BoundsError(ws.fixed_capital,chosen_index))
    selected_K=ws.fixed_capital[chosen_index]
    selected_X=ws.fixed_investment[chosen_index]
    selected_X>0.0 ||
        throw(DomainError(selected_X,
            "nonpositive X cannot enter original log-investment closure"))
    sumX=0.0
    investing=0
    divesting=0
    @inbounds for i in 1:n
        d=ws.depreciated[i]
        lambda_i=state.firms.lambda[i]
        if selected_K<=ws.switches[2i-1]
            newk=d/lambda_i
            divesting+=1
        elseif selected_K<=ws.switches[2i]
            newk=d
        else
            newk=d*lambda_i
            investing+=1
        end
        ws.next_capital[i]=newk
        sumX+=newk-d
    end
    check_K=aggregate_capital(
        state.firms.next_productivity,ws.next_capital,p.alpha,rho
    )
    check_X=sumX/n
    abs(check_K-selected_K)<=1e-10*max(1.0,selected_K) ||
        error("Diagnostic selected capital fails independent CES check")
    abs(check_X-selected_X)<=1e-10*max(1.0,abs(selected_X)) ||
        error("Diagnostic selected investment fails independent accounting")
    old_capital=state.firms.capital
    state.firms.capital=ws.next_capital
    ws.next_capital=old_capital
    advance_productivity!(state.firms)
    state.K=selected_K
    state.mu_prior=cal.steady_state.mu*
        exp(state.re.h*log(selected_K/cal.steady_state.K))
    state.period+=1
    return Dict{String,Any}(
        "K"=>selected_K,"X"=>selected_X,
        "epsilon"=>log(selected_X/expected_X),
        "investing_firms"=>investing,"divesting_firms"=>divesting,
        "capital_check_residual"=>check_K-selected_K,
        "investment_check_residual"=>check_X-selected_X,
    )
end

function _diagnostic_selection_index(ranks, rule::String, time::Int, deviation::Int)
    if rule=="published_investment_ES1"
        return ranks.invest_best
    elseif rule=="capital_distance"
        return ranks.capital_best
    elseif rule=="maximum_positive_capital"
        return ranks.max_positive
    elseif rule=="second_best_at_one_date"
        return time==deviation ? ranks.second_best : ranks.invest_best
    else
        error("Unsupported diagnostic selection rule: "*rule)
    end
end

function _selection_path_trial(
    rule::String;deviation::Int=-1,horizon::Int=20,save_full::Bool=false,
    baseline_check::Bool=false
)
    state=initialize_simulation(;N=110000,base_seed=20261101)
    check_state=baseline_check ?
        initialize_simulation(;N=110000,base_seed=20261101) : nothing
    periods=Vector{Dict{String,Any}}()
    completed=0
    failure="none"
    disagreement_dates=Int[]
    for t in 0:horizon-1
        rec,ranks,ex,ek=_period_selection_audit(state)
        if ranks.capital_best!=ranks.invest_best
            push!(disagreement_dates,t)
        end
        chosen=_diagnostic_selection_index(ranks,rule,t,deviation)
        rec["chosen_index"]=chosen
        rec["rule"]=rule
        rec["one_time_deviation_date"]=deviation
        if chosen==0
            failure="No positive-X fixed point available for selection"
            rec["failure"]=failure
            push!(periods,rec)
            break
        end
        if state.workspace.fixed_investment[chosen]<=0.0
            failure="Chosen capital-based equilibrium has nonpositive X"
            rec["failure"]=failure
            push!(periods,rec)
            break
        end
        before_k=baseline_check ? copy(state.firms.capital) : nothing
        advanced=_advance_with_diagnostic_selection!(state,chosen,ex)
        rec["selected_state_transition"]=advanced
        if baseline_check
            # Independently execute full production GE+ES1 and require that
            # the diagnostic copy follows the same state at every date.
            selected_production=simulate_period!(check_state)
            ranks.invest_best==chosen ||
                error("Benchmark comparison called on non-ES1 rule")
            selected_production.K_next_firm==advanced["K"] ||
                error("Diagnostic/production K mismatch")
            isapprox(selected_production.X,advanced["X"];atol=1e-12,rtol=1e-12) ||
                error("Diagnostic/production X mismatch")
            state.firms.capital==check_state.firms.capital ||
                error("Diagnostic/production firm capital mismatch")
            state.firms.productivity==check_state.firms.productivity ||
                error("Diagnostic/production productivity mismatch")
            state.K==check_state.K ||
                error("Diagnostic/production GE state K mismatch")
            isapprox(state.mu_prior,check_state.mu_prior;atol=1e-10,rtol=1e-12) ||
                error("Diagnostic/production stable continuation mismatch")
            rec["baseline_production_transition_exactly_verified"]=true
        end
        push!(periods,rec)
        completed+=1
    end
    return Dict{String,Any}(
        "rule"=>rule,"deviation"=>deviation,
        "requested_periods"=>horizon,"completed_periods"=>completed,
        "failure"=>failure,
        "capital_and_investment_selection_disagreement_dates"=>
            disagreement_dates,
        "periods"=>periods,
        "completed"=>completed==horizon,
        "diagnostic_not_published_ES1"=>rule!="published_investment_ES1",
    )
end

function equilibrium_path_dependence_audit()
    # The full baseline path is accompanied by independent checks
    # against the UNMODIFIED simulate_period! production implementation.
    baseline=_selection_path_trial(
        "published_investment_ES1";horizon=20,
        baseline_check=true
    )
    first=baseline["periods"][1]
    first["investment_ES1_index"]==first["total_fixed_points"] ||
        error("Known benchmark t0 ES1 must select maximum-capital point")
    isapprox(first["investment_ES1_K"],8.546731346350889;
             atol=1e-10,rtol=1e-12) ||
        error("Benchmark ES1 initial capital regression changed")
    isapprox(first["investment_ES1_X"],0.06936827009543808;
             atol=1e-10,rtol=1e-12) ||
        error("Benchmark ES1 initial investment regression changed")
    trials=Dict{String,Any}(
        "capital_distance"=>_selection_path_trial("capital_distance"),
        "maximum_positive_capital"=>_selection_path_trial(
            "maximum_positive_capital"
        ),
    )
    for tstar in 0:4
        label="second_best_one_time_t"*string(tstar)
        trials[label]=_selection_path_trial(
            "second_best_at_one_date";deviation=tstar
        )
    end
    # Verify identical prefixes of all realized macro paths before the
    # first deliberate deviation, and identical evolution if selection
    # coincides with benchmark at that date.
    for tstar in 0:4
        trial=trials["second_best_one_time_t"*string(tstar)]
        for t in 0:(tstar-1)
            reference=baseline["periods"][t+1]
            candidate=trial["periods"][t+1]
            reference["chosen_index"]==candidate["chosen_index"] ||
                error("Pre-deviation equilibrium selection differs")
            reftr=reference["selected_state_transition"]
            cftr=candidate["selected_state_transition"]
            reftr["K"]==cftr["K"] && reftr["X"]==cftr["X"] ||
                error("Pre-deviation capital or investment differs")
        end
        deviated=trial["periods"][tstar+1]
        original_at_deviation=baseline["periods"][tstar+1]
        deviated["chosen_index"]!=original_at_deviation["chosen_index"] ||
            error("Second-best experiment did not actually deviate")
        devtr=deviated["selected_state_transition"]
        basetr=original_at_deviation["selected_state_transition"]
        devtr["K"]!=basetr["K"] ||
            error("Second-best experiment did not change aggregate capital")
    end
    return Dict{String,Any}(
        "interpretation"=>"Purely diagnostic alternate fixed-point selection, original iid shocks, original {-1,0,1} and Appendix A macro recursion. Not a new RE equilibrium.",
        "baseline"=>baseline,
        "alternative_paths"=>trials,
        "baseline_t0_highest_capital_regression_passed"=>true,
        "baseline_diagnostic_transition_verified_against_original_src"=>true,
        "common_future_productivity_rng_seed"=>20261101,
    )
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
        four_way = four_way_productivity_effects(
            state.firms, cal, Ke, no_new_productivity, at_expected
        )
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
            "four_way_productivity_at_expected_K"=>four_way,
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
        if t==0
            # Exact-in-gap moments and Gaussian quadrature for future iid
            # productivity. Do not imply that finite-N forecast K is fixed.
            gapmom_24 = initial_uniform_gap_expected_moments(
                state.firms, cal, Ke; nodes=24
            )
            gapmom_36 = initial_uniform_gap_expected_moments(
                state.firms, cal, Ke; nodes=36
            )
            info["initial_uniform_gap_24_nodes"] = gapmom_24
            info["initial_uniform_gap_36_nodes"] = gapmom_36
            info["initial_uniform_gap_convergence_X"] =
                gapmom_36["expected_net_investment"]-
                gapmom_24["expected_net_investment"]
            info["initial_uniform_gap_expected_X_minus_continuum_forecast"] =
                gapmom_36["expected_net_investment"]-expected.X_next
            continuum_low=initial_continuum_consistency_audit(
                cal, Ke;lambda_nodes=32,shock_nodes=40
            )
            continuum_high=initial_continuum_consistency_audit(
                cal, Ke;lambda_nodes=48,shock_nodes=64
            )
            info["unconditional_initial_continuum_32x40"]=continuum_low
            info["unconditional_initial_continuum_48x64"]=continuum_high
            gaussian_exact_48=exact_initial_continuum_consistency_audit(
                cal,Ke;lambda_nodes=48
            )
            gaussian_exact_64=exact_initial_continuum_consistency_audit(
                cal,Ke;lambda_nodes=64
            )
            info["exact_normal_continuum_48_lambda_nodes"]=gaussian_exact_48
            info["exact_normal_continuum_64_lambda_nodes"]=gaussian_exact_64
            info["exact_normal_continuum_lambda_convergence_Phi"]=
                gaussian_exact_64["Phi_at_Kbar_minus_Kbar"]-
                gaussian_exact_48["Phi_at_Kbar_minus_Kbar"]
            info["continuum_quadrature_convergence_Phi_at_Kbar"]=
                continuum_high["Phi_at_Kbar"]-
                continuum_low["Phi_at_Kbar"]
            info["continuum_quadrature_convergence_X_at_Kbar"]=
                continuum_high["X_at_Kbar"]-
                continuum_low["X_at_Kbar"]
        end
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
    output["support_condition_counterfactual"] =
        support_consistent_comparison(;horizon=20)
    output["equilibrium_selection_path_dependence"] =
        equilibrium_path_dependence_audit()
    output["preconditioned_distribution_counterfactual"] =
        preconditioned_distribution_comparison()
    output["exogenous_fixed_K_firm_distribution"] =
        fixed_K_distribution_diagnostic(
            cal;n=n,seed=20261101,horizon=256
        )
    path = normpath(joinpath(@__DIR__,"..","results","initial_collapse_audit.toml"))
    mkpath(dirname(path))
    open(path,"w") do io
        TOML.print(io,output)
    end
    println("Audit saved: ",path)
    println("One-sided-condition violations: ",violated," of ",n)
    lambda_study=output["support_condition_counterfactual"]
    lambda_A=lambda_study["case_A_initial_equilibrium"]
    lambda_B=lambda_study["case_B_initial_equilibrium"]
    lambda_path=lambda_study["case_B_simulation"]
    precond=output["preconditioned_distribution_counterfactual"]
    pa=precond["first_period"]["caseA"]
    pb=precond["first_period"]["caseB"]
    pc=precond["caseB_twenty_periods"]
    println("Preconditioned initial A: K0=",pa["K0"],
        "; expected X=",pa["expected_X"],
        "; fixed points=",pa["fixed_points"],
        "; ES1 X=",get(pa,"ES1_X","none"))
    println("Preconditioned initial B: K0=",pb["K0"],
        "; expected X=",pb["expected_X"],
        "; reaction residual=",pb["reaction_residual_at_expected_K"],
        "; fixed points=",pb["fixed_points"],
        "; positive X fixed points=",pb["positive_investment_fixed_points"],
        "; ES1 X=",get(pb,"ES1_X","none"),
        "; ES1 K=",get(pb,"ES1_K","none"))
    println("Preconditioned B gap outside [0,1): ",
        precond["first_period"]["caseB_preconditioned_gap"]["fraction_gap_outside_0_1"])
    println("Preconditioned B GE completed ",
        pc["completed_periods"],"/",pc["requested_periods"],
        " periods; failure=",pc["failure"])
    ranking=output["equilibrium_selection_path_dependence"]
    base_path=ranking["baseline"]
    println("Selection path baseline: completed=",base_path["completed_periods"],
        "/20; failure=",base_path["failure"],
        "; capital-vs-investment disagreement dates=",
        base_path["capital_and_investment_selection_disagreement_dates"])
    for row in base_path["periods"]
        println("Selection baseline t=",row["date"],
            "; fixed points=",row["total_fixed_points"],
            "; positive X=",row["positive_investment_fixed_points"],
            "; selected capital rank=",
            get(row,"investment_ES1_capital_rank_ascending","none"),
            "/",row["total_fixed_points"],
            "; capital-distance index=",row["capital_distance_index"],
            "; X=",get(row,"investment_ES1_X","none"))
    end
    for label in sort(collect(keys(ranking["alternative_paths"])))
        trial=ranking["alternative_paths"][label]
        start_t=trial["deviation"]
        changed_info=if start_t>=0
            base_t=base_path["periods"][start_t+1]
            alt_t=trial["periods"][start_t+1]
            bt=base_t["selected_state_transition"]
            at=alt_t["selected_state_transition"]
            string("; deviation date=",start_t,
                "; delta selected K=",at["K"]-bt["K"],
                "; delta selected X=",at["X"]-bt["X"])
        else
            ""
        end
        final_candidate=trial["periods"][end]
        println("Selection alternative ",label,": completed=",
            trial["completed_periods"],"/20; failure=",trial["failure"],
            "; disagreements=",
            trial["capital_and_investment_selection_disagreement_dates"],
            changed_info,
            "; final date=",final_candidate["date"],
            "; final fixed points=",final_candidate["total_fixed_points"],
            "; final positive-X points=",
            final_candidate["positive_investment_fixed_points"],
            "; final max X=",final_candidate["maximum_X"])
    end
    println("Support-only Case A: t0 ES1 X=",lambda_A["ES1_X"],
        "; K=",lambda_A["ES1_K"],
        "; first-period residual=",lambda_A["reaction_residual_at_expected_K"])
    println("Support-only Case B: changed firms=",
        lambda_study["case_B_distribution"]["modified_firms"],
        "; support violations=",lambda_study["case_B_distribution"]["violations"],
        "; t0 ES1 X=",get(lambda_B,"ES1_X","none"),
        "; K=",get(lambda_B,"ES1_K","none"),
        "; first-period residual=",lambda_B["reaction_residual_at_expected_K"])
    println("Support-only Case B: completed periods=",
        lambda_path["completed_periods"],"/",lambda_path["requested_periods"],
        "; failure=",lambda_path["failure"])
    sample_A=lambda_study["case_A_distribution"]
    sample_B=lambda_study["case_B_distribution"]
    direct_lambda=lambda_study["direct_modified_firm_effect_at_A_expected_K"]
    println("Support-only lambda distribution A: mean=",sample_A["mean"],
        "; sd=",sample_A["sd"],
        "; min=",sample_A["minimum"],
        "; max=",sample_A["maximum"])
    println("Support-only lambda distribution B: mean=",sample_B["mean"],
        "; sd=",sample_B["sd"],
        "; min=",sample_B["minimum"],
        "; max=",sample_B["maximum"])
    println("Support-only initial A: K0=",lambda_A["K0"],
        "; expected X=",lambda_A["expected_X"],
        "; expected K=",lambda_A["expected_K"],
        "; fixed points=",lambda_A["fixed_points"],
        "; positive-X fixed points=",lambda_A["positive_investment_fixed_points"])
    println("Support-only initial B: K0=",lambda_B["K0"],
        "; expected X=",lambda_B["expected_X"],
        "; expected K=",lambda_B["expected_K"],
        "; fixed points=",lambda_B["fixed_points"],
        "; positive-X fixed points=",lambda_B["positive_investment_fixed_points"])
    println("Support-only changed 673 firms at fixed K_A: direct X effect=",
        direct_lambda["X_direct_change_at_fixed_K"],
        "; direct CES power effect=",
        direct_lambda["CES_power_direct_change_at_fixed_K"],
        "; direct aggregate K effect=",
        direct_lambda["Phi_direct_change_at_fixed_K"])
    cross=output["exogenous_fixed_K_firm_distribution"]
    first_state=cross["snapshots"][1]
    final_state=cross["snapshots"][end]
    println("Fixed-K three-choice distribution: initial vs date 256 ",
        "unwrapped-gap-outside=",
        first_state["fraction_actual_gap_outside_0_1"]," vs ",
        final_state["fraction_actual_gap_outside_0_1"],
        "; fractional-gap max 10-bin deviation=",
        first_state["max_fractional_gap_bin_error_from_uniform"]," vs ",
        final_state["max_fractional_gap_bin_error_from_uniform"])
    println("Fixed-K final CES-minus-Kbar=",
        final_state["CES_minus_exogenous_Kbar"],
        "; two recent block mean X=",cross["last_two_64_period_mean_gross_investment"],
        "; two recent block mean capital drift=",
        cross["last_two_64_period_mean_change_in_firm_capital"])
    for rec in periods
        if rec["t"]==0
            four = rec["four_way_productivity_at_expected_K"]
            moments = rec["initial_uniform_gap_36_nodes"]
            println("t=0 four-way reaction residuals A/B/C/D: ",
                four["residual_A"], " / ",four["residual_B"], " / ",
                four["residual_C"], " / ",four["residual_D"])
            println("t=0 Shapley threshold effect: ",
                four["threshold_effect_shapley"],
                "; direct CES effect: ",four["CES_weight_effect_shapley"])
            println("t=0 uniform-gap and iid-productivity expected X (fixed K): ",
                moments["expected_net_investment"],
                "; 36 minus 24 node error: ",
                rec["initial_uniform_gap_convergence_X"])
            continuum=rec["unconditional_initial_continuum_48x64"]
            println("t=0 unconditional initial-continuum Phi(Kbar)-Kbar: ",
                continuum["Phi_at_Kbar_minus_Kbar"],
                "; X(Kbar)-Xbar: ",continuum["X_at_Kbar_minus_Xbar"],
                "; roots in searched range: ",
                continuum["root_capitals_in_tested_intervals"],
                "; investment at roots: ",continuum["root_investments"])
            println("t=0 unconditional continuum quadrature convergence Phi: ",
                rec["continuum_quadrature_convergence_Phi_at_Kbar"],
                "; X: ",rec["continuum_quadrature_convergence_X_at_Kbar"])
            exact=rec["exact_normal_continuum_64_lambda_nodes"]
            println("t=0 exact-Gaussian continuum Phi(Kbar)-Kbar: ",
                exact["Phi_at_Kbar_minus_Kbar"],
                "; X(Kbar)-Xbar: ",exact["X_at_Kbar_minus_Xbar"],
                "; roots: ",exact["root_capitals_in_tested_intervals"],
                "; root X: ",exact["root_investments"])
            println("t=0 Gaussian-analytic lambda quadrature convergence Phi: ",
                rec["exact_normal_continuum_lambda_convergence_Phi"])
        end
        println("t=",rec["t"]," expected K=",rec["expected_K"],
            " Phi(expected K)-expected K=",rec["reaction_residual_at_expected_K"],
            " selected X=",get(rec,"X_selected","no positive-X fixed point"),
            " fixed points=",rec["n_fixed"])
    end
    return nothing
end

main()
