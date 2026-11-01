#!/usr/bin/env julia
# Book-motivated two-choice counterfactual for Nirei (2015), Issue #5.
# This is NOT the published 2015 three-choice benchmark.
# From Nirei-2015/:
#   julia --project=. scripts/run_two_choice.jl --smoke
#   julia --project=. scripts/run_two_choice.jl 6 110000 20261101 1 hybrid
#   julia --project=. scripts/run_two_choice.jl 6 110000 20261101 1 linear
#   julia --project=. scripts/run_two_choice.jl 20 110000 20261101
#   julia --project=. scripts/run_two_choice.jl 1100 110000 20261101

using Nirei2015
using Statistics: mean, std, cor
using TOML

# Only investment or inaction is permitted:
# k_{i,t+1} = (1-delta)*k_{i,t} * lambda_i^kappa, kappa in {0,1}.
# The investment threshold is the original upper switching boundary d_i/b_i.
# At equality the firm remains inactive, exactly as in the original code.
# No divestment, clipping, shock rejection, or additional firm adjustment.
function enumerate_two_choice(firms::FirmState, calibration::CalibrationResult)
    firms.next_ready ||
        throw(ArgumentError("draw_productivity! must precede two-choice enumeration"))
    p = calibration.parameters
    rho = calibration.constants.rho
    n = length(firms.capital)
    n > 0 || throw(ArgumentError("N must be positive"))
    thresholds = Vector{Float64}(undef, n)
    depreciated = Vector{Float64}(undef, n)
    base_power = Vector{Float64}(undef, n)
    power_sum = 0.0
    for i in 1:n
        ai = firms.next_productivity[i]
        lambda = firms.lambda[i]
        d = (1.0 - p.delta) * firms.capital[i]
        b = threshold_coefficient(lambda, ai, calibration)
        threshold = d / b
        base = ai^(rho / p.alpha) * d^rho / n
        0.0 < threshold < Inf ||
            throw(DomainError(threshold, "invalid two-choice investment threshold"))
        isfinite(base) && base > 0.0 ||
            throw(DomainError(base, "invalid CES weight"))
        thresholds[i] = threshold
        depreciated[i] = d
        base_power[i] = base
        power_sum += base
    end

    perm = sortperm(thresholds)
    capital_roots = Float64[]
    investment_roots = Float64[]
    invrho = inv(rho)
    investment_sum = 0.0
    left = 0.0
    j = 1
    while j <= n
        boundary = thresholds[perm[j]]
        candidate = power_sum^invrho
        if left < candidate <= boundary
            push!(capital_roots, candidate)
            push!(investment_roots, investment_sum)
        end
        while j <= n && thresholds[perm[j]] == boundary
            i = perm[j]
            lambda = firms.lambda[i]
            power_sum += base_power[i] * (lambda^rho - 1.0)
            investment_sum += depreciated[i] * (lambda - 1.0) / n
            j += 1
        end
        left = boundary
    end
    candidate = power_sum^invrho
    if left < candidate
        push!(capital_roots, candidate)
        push!(investment_roots, investment_sum)
    end
    isempty(capital_roots) &&
        throw(ErrorException("No two-choice aggregate capital fixed point"))
    return (
        K=capital_roots, X=investment_roots,
        thresholds=thresholds, depreciated=depreciated,
    )
end

# Preserve the original Section 3.4 log-investment ES1 over X>0.
# Two-choice prevents negative X but need not guarantee a strictly positive X.
function select_two_choice_es1(roots, expected_X::Float64, expected_K::Float64)
    expected_X > 0.0 && expected_K > 0.0 ||
        throw(DomainError((expected_X, expected_K), "positive forecasts required"))
    best = 0
    best_score = Inf
    best_distance = Inf
    positive_count = 0
    for j in eachindex(roots.K)
        xj = roots.X[j]
        if xj > 0.0
            positive_count += 1
            score = abs(log(xj / expected_X))
            distance = abs(log(roots.K[j] / expected_K))
            if score < best_score ||
               (score == best_score &&
                (distance < best_distance ||
                 (distance == best_distance &&
                  (best == 0 || roots.K[j] < roots.K[best]))))
                best = j
                best_score = score
                best_distance = distance
            end
        end
    end
    positive_count > 0 ||
        throw(DomainError((
            fixed_points=length(roots.K),
            minimum_X=minimum(roots.X),
            maximum_X=maximum(roots.X),
            expected_X=expected_X,
        ), "Two-choice equilibrium has no strictly positive X; log ES1 undefined"))
    return (index=best, score=best_score, n_positive=positive_count)
end

# Appendix A first-order realized macro allocation, NOT an exact nonlinear GE solve.
# Equations (27),(28),(30),(33) imply the expressions below.
# The precommitted interest-rate proxy is the *conditional forecast* of R;
# it must not be interpreted as an ex-post household Euler identity.
function linear_realized_macroeconomy(
    state::SimulationState, K_prior::Float64, X::Float64,
    logmu::Float64, forecast::ExpectedEquilibrium,
)
    p = state.calibration.parameters
    c = state.calibration.constants
    s = state.calibration.steady_state
    coeff = state.re.coefficients
    logk = log(K_prior / s.K)
    logx = log(X / s.X)
    logl = (s.a_x * logx - p.alpha * logk -
            (s.a_c / p.sigma) * logmu) / coeff.A1
    logc = (-logmu + (p.sigma - 1.0) * (c.c_L / s.a_c) * logl) / p.sigma
    logy = p.alpha * logk + (1.0 - p.alpha) * logl
    logw = logc + (p.zeta - 1.0 + c.c_L / s.a_c) * logl
    all(isfinite, (logk, logx, logl, logc, logy, logw, logmu)) ||
        throw(DomainError((logk, logx, logl, logc, logy, logw, logmu),
                          "nonfinite Appendix A log-deviations"))
    residuals = (
        logy - p.alpha * logk - (1.0 - p.alpha) * logl,
        logy - s.a_c * logc - s.a_x * logx,
        logmu + p.sigma * logc -
            (p.sigma - 1.0) * (c.c_L / s.a_c) * logl,
        logl - coeff.eta_L * (logw + logmu / p.sigma),
        logw - logc - (p.zeta - 1.0 + c.c_L / s.a_c) * logl,
    )
    linear_residual = maximum(abs, residuals)
    linear_residual < 1e-10 ||
        throw(ErrorException("Appendix A equation residual = $linear_residual"))
    Y = s.Y * exp(logy)
    C = s.C * exp(logc)
    L = s.L * exp(logl)
    w = s.w * exp(logw)
    R = forecast.R_next
    all(z -> isfinite(z) && z > 0.0, (Y, C, L, w, R)) ||
        throw(DomainError((Y, C, L, w, R),
                          "invalid level reconstructed from first-order approximation"))
    return (
        Y=Y, C=C, L=L, w=w, R=R,
        resource_residual=(Y - C - X)/Y,
        linear_residual=linear_residual,
    )
end

# For the hybrid branch the macro model is solved at nonlinear levels.
# Evaluate the *first-order* approximation residual as a diagnostic only.
function first_order_macro_residual(
    state::SimulationState, K_prior::Float64, X::Float64,
    logmu::Float64, allocation,
)
    p = state.calibration.parameters
    c = state.calibration.constants
    s = state.calibration.steady_state
    coeff = state.re.coefficients
    logk = log(K_prior / s.K)
    logx = log(X / s.X)
    logl = log(allocation.L / s.L)
    logc = log(allocation.C / s.C)
    logy = log(allocation.Y / s.Y)
    logw = log(allocation.w / s.w)
    return maximum(abs, (
        logy - p.alpha * logk - (1.0-p.alpha) * logl,
        logy - s.a_c*logc - s.a_x*logx,
        logmu + p.sigma*logc - (p.sigma-1.0)*(c.c_L/s.a_c)*logl,
        logl - coeff.eta_L*(logw + logmu/p.sigma),
        logw - logc - (p.zeta-1.0+c.c_L/s.a_c)*logl,
    ))
end

# Keep the original firm shocks, two-choice investment rule and ES1.
# Only the macroeconomic continuation differs across variants; the hybrid
# preserves the prior exact-level solver and the linear variant uses Appendix A.
function simulate_two_choice_period!(state::SimulationState, continuation::Symbol=:hybrid)
    continuation in (:hybrid, :linear) ||
        throw(ArgumentError("unknown macro continuation: $continuation"))
    p = state.calibration.parameters
    s = state.calibration.steady_state
    kp = state.K
    prior_mu = state.mu_prior
    t = state.period
    try
        forecast = compute_expected_equilibrium(
            state.re, state.calibration, kp, prior_mu
        )
        draw_productivity!(state.firms, state.rngs, state.calibration)
        roots = enumerate_two_choice(state.firms, state.calibration)
        expected_K = (1.0 - p.delta) * kp +
                     state.calibration.threshold.A^(state.calibration.constants.rho - 1.0) /
                     state.calibration.constants.rho * forecast.X_next
        es1 = select_two_choice_es1(roots, forecast.X_next, expected_K)
        chosen_K = roots.K[es1.index]
        chosen_X = roots.X[es1.index]
        n = length(state.firms.capital)
        next_capital = state.workspace.next_capital
        investing = 0
        investment_check = 0.0
        @inbounds for i in 1:n
            d = roots.depreciated[i]
            k = if chosen_K > roots.thresholds[i]
                investing += 1
                d * state.firms.lambda[i]
            else
                d
            end
            next_capital[i] = k
            investment_check += k - d
        end
        K_check = aggregate_capital(
            state.firms.next_productivity, next_capital,
            p.alpha, state.calibration.constants.rho
        )
        X_check = investment_check / n
        K_residual = K_check - chosen_K
        X_residual = X_check - chosen_X
        abs(K_residual) <= 1e-10 * max(1.0, abs(chosen_K)) ||
            throw(ErrorException("two-choice CES fixed-point residual = $K_residual"))
        abs(X_residual) <= 1e-10 * max(1.0, abs(chosen_X)) ||
            throw(ErrorException("two-choice investment residual = $X_residual"))

        epsilon = log(chosen_X / forecast.X_next)
        logmu_prior = log(prior_mu / s.mu)
        logK_next = log(chosen_K / s.K)
        # Hybrid reproduces the existing numerical continuation.
        # Linear explicitly reanchors the *realized* jump variable to the
        # local stable manifold after exact finite-N CES aggregation.
        # The latter is a model approximation, NOT an exact global RE result.
        logmu = continuation === :linear ?
            state.re.h * logK_next :
            forecast.logmu_next + state.re.h * p.delta * epsilon
        mu = s.mu * exp(logmu)
        isfinite(mu) && mu > 0.0 ||
            throw(DomainError(mu, "invalid two-choice continuation marginal utility"))
        realized_macro = if continuation === :linear
            linear_realized_macroeconomy(state, kp, chosen_X, logmu, forecast)
        else
            Nirei2015._realized_macroeconomy(
                state.calibration, kp, chosen_X, logmu, logmu_prior
            )
        end
        macro_linear_residual = continuation === :linear ?
            realized_macro.linear_residual :
            first_order_macro_residual(state, kp, chosen_X, logmu, realized_macro)
        stable_gap = logmu - state.re.h * logK_next
        investment_capital_gap = logK_next - forecast.logK_next -
                                 p.delta * epsilon
        all(isfinite, (macro_linear_residual, stable_gap,
                       investment_capital_gap)) ||
            throw(ErrorException("nonfinite macro continuation diagnostics"))
        record = (
            t=t, K_prior=kp, K_next=chosen_K, expected_K=expected_K,
            X=chosen_X, expected_X=forecast.X_next, epsilon=epsilon,
            Y=realized_macro.Y, C=realized_macro.C, L=realized_macro.L,
            w=realized_macro.w, R=realized_macro.R, mu=mu,
            logK=log(chosen_K), logX=log(chosen_X),
            logY=log(realized_macro.Y), logC=log(realized_macro.C), logL=log(realized_macro.L),
            n_fixed=length(roots.K), n_positive=es1.n_positive,
            n_investing=investing, n_inactive=n-investing, n_divesting=0,
            K_residual=K_residual, X_residual=X_residual,
            resource_residual=realized_macro.resource_residual,
            macro_linear_residual=macro_linear_residual,
            logmu=logmu, stable_gap=stable_gap,
            expected_logK=forecast.logK_next,
            expected_logmu=forecast.logmu_next,
            capital_linear_gap=investment_capital_gap,
        )
        # Commit firm capital and productivity only after all checks pass.
        old_capital = state.firms.capital
        state.firms.capital = next_capital
        state.workspace.next_capital = old_capital
        advance_productivity!(state.firms)
        state.K = chosen_K
        state.mu_prior = mu
        state.period += 1
        return record
    catch err
        err isa InterruptException && rethrow()
        throw(SimulationFailure(t, kp, prior_mu, sprint(showerror, err)))
    end
end

function smoke_test()
    state = initialize_simulation(; N=8, base_seed=20261101)
    draw_productivity!(state.firms, state.rngs, state.calibration)
    roots = enumerate_two_choice(state.firms, state.calibration)
    p = state.calibration.parameters
    rho = state.calibration.constants.rho
    firms = state.firms
    n = length(firms.capital)
    reference = Tuple{Float64,Float64}[]
    for mask in 0:((1 << n)-1)
        capitals = Vector{Float64}(undef, n)
        investment = 0.0
        for i in 1:n
            d = (1.0-p.delta)*firms.capital[i]
            invest = ((mask >> (i-1)) & 1) == 1
            capitals[i] = d * (invest ? firms.lambda[i] : 1.0)
            investment += capitals[i] - d
        end
        K = aggregate_capital(firms.next_productivity, capitals, p.alpha, rho)
        is_equilibrium = all(
            (((mask >> (i-1)) & 1) == 1) ==
            (K > roots.thresholds[i]) for i in 1:n
        )
        is_equilibrium && push!(reference, (K, investment/n))
    end
    sort!(reference; by=first)
    length(reference) == length(roots.K) ||
        error("independent 2^N oracle: wrong root count")
    for j in eachindex(reference)
        isapprox(reference[j][1], roots.K[j]; rtol=1e-11, atol=1e-11) ||
            error("independent oracle: capital root mismatch")
        isapprox(reference[j][2], roots.X[j]; rtol=1e-11, atol=1e-11) ||
            error("independent oracle: investment mismatch")
    end
    all(x -> x >= 0.0, roots.X) ||
        error("two-choice net investment cannot be negative")
    println("PASS: two-choice N=8 exhaustive 2^N policy oracle; roots=", length(reference))
    # Nontrivial first-order resource, production, utility and labor FOC
    # regression, separate from the exhaustive firm-policy oracle above.
    fresh = initialize_simulation(; N=8, base_seed=20261101)
    steady = fresh.calibration.steady_state
    expected = compute_expected_equilibrium(
        fresh.re, fresh.calibration, steady.K, steady.mu
    )
    K_test = 1.01 * steady.K
    X_test = 1.05 * steady.X
    logmu_test = fresh.re.h * log(0.99)
    allocation = linear_realized_macroeconomy(
        fresh, K_test, X_test, logmu_test, expected
    )
    allocation.linear_residual < 1e-12 ||
        error("Appendix A first-order equations failed at nonsteady test state")
    println("PASS: Appendix A linear equations on nonsteady test state")

end

function write_csv(path::String, records)
    isempty(records) && return
    names = propertynames(first(records))
    open(path, "w") do io
        println(io, join(string.(names), ","))
        for rec in records
            println(io, join((string(getproperty(rec, name)) for name in names), ","))
        end
    end
end

function main()
    if length(ARGS) >= 1 && ARGS[1] == "--smoke"
        smoke_test()
        return
    end
    T = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 6
    N = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 110000
    seed = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 20261101
    replication = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 1
    continuation = length(ARGS) >= 5 ? Symbol(ARGS[5]) : :hybrid
    continuation in (:hybrid, :linear) ||
        throw(ArgumentError("continuation must be hybrid or linear"))
    T > 0 && N > 0 && replication > 0 ||
        throw(ArgumentError("T, N and replication must be positive"))
    state = initialize_simulation(; N=N, base_seed=seed, replication=replication)
    records = NamedTuple[]
    failure = nothing
    for _ in 1:T
        try
            push!(records, simulate_two_choice_period!(state, continuation))
        catch err
            failure = sprint(showerror, err)
            break
        end
    end
    root = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(root)
    stem = "two_choice_rep" * lpad(string(replication), 2, '0') * "_" * string(continuation)
    csv_path = joinpath(root, stem * ".csv")
    summary_path = joinpath(root, stem * "_summary.toml")
    write_csv(csv_path, records)
    summary = Dict{String,Any}(
        "status" => Dict(
            "completed" => failure === nothing && length(records) == T,
            "completed_periods" => length(records),
            "requested_periods" => T,
            "failure" => something(failure, "none"),
        ),
        "configuration" => Dict(
            "model" => "book-motivated two-choice kappa in {0,1}",
            "macro_variant" => string(continuation),
            "not_published_2015_replication" => true,
            "firms" => N, "seed" => seed, "replication" => replication,
            "selection" => "original Section 3.4 ES1 over X>0",
            "continuation" => continuation === :linear ?
                "realized K reanchored to local saddle; log-linear Appendix A variables" :
                "existing first-order mu shock plus exact nonlinear household levels",
            "burnin" => 100,
        ),
    )
    if !isempty(records)
        r = first(records)
        summary["first_period"] = Dict(
            "K_next" => r.K_next, "X" => r.X,
            "expected_X" => r.expected_X, "n_fixed" => r.n_fixed,
            "n_positive" => r.n_positive, "n_investing" => r.n_investing,
            "capital_residual" => r.K_residual,
            "stable_gap" => r.stable_gap,
            "max_linear_macro_residual" => r.macro_linear_residual,
        )
        r = last(records)
        summary["last_period"] = Dict(
            "K_next" => r.K_next, "X" => r.X, "n_fixed" => r.n_fixed,
            "mu" => r.mu, "stable_gap" => r.stable_gap,
        )
    end
    if !isempty(records)
        summary["diagnostics"] = Dict(
            "maximum_abs_stable_gap" => maximum(r -> abs(r.stable_gap), records),
            "maximum_macro_linear_residual" => maximum(r -> r.macro_linear_residual, records),
            "maximum_abs_capital_linear_gap" => maximum(r -> abs(r.capital_linear_gap), records),
            "last_capital" => last(records).K_next,
        )
    end
    # Do not calculate post-burn-in benchmark moments from truncated failures.
    if failure === nothing && length(records) > 101
        keep = records[101:end]
        summary["moments"] = Dict(
            "sd_logY_pct" => 100*std(r.logY for r in keep),
            "sd_logC_pct" => 100*std(r.logC for r in keep),
            "sd_logX_pct" => 100*std(r.logX for r in keep),
            "sd_logL_pct" => 100*std(r.logL for r in keep),
            "sd_logK_pct" => 100*std(r.logK for r in keep),
            "corr_Y_C" => cor([r.logY for r in keep], [r.logC for r in keep]),
            "corr_Y_X" => cor([r.logY for r in keep], [r.logX for r in keep]),
            "corr_Y_L" => cor([r.logY for r in keep], [r.logL for r in keep]),
            "corr_Y_K" => cor([r.logY for r in keep], [r.logK for r in keep]),
        )
    end
    open(summary_path, "w") do io
        TOML.print(io, summary)
    end
    println("Two-choice counterfactual: T=$T N=$N seed=$seed replication=$replication continuation=$continuation")
    println("Completed ", length(records), "/", T, " periods")
    if !isempty(records)
        r = first(records)
        println("t=0: X=", r.X, " Knext=", r.K_next,
                " fixed points=", r.n_fixed, " positive-X=", r.n_positive,
                " investing=", r.n_investing)
        r = last(records)
        println("last: t=", r.t, " X=", r.X, " Knext=", r.K_next,
                " mu=", r.mu, " stable_gap=", r.stable_gap)
    end
    println("CSV: ", csv_path)
    println("Summary: ", summary_path)
    if failure !== nothing
        error("Two-choice simulation failed: " * failure)
    end
    if haskey(summary, "moments")
        println("Post-burn-in moments: ", summary["moments"])
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
