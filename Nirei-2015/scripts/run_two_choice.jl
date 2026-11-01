#!/usr/bin/env julia
# Book-motivated two-choice counterfactual for Nirei (2015), Issue #5.
# This is NOT the published 2015 three-choice benchmark.
# From Nirei-2015/:
#   julia --project=. scripts/run_two_choice.jl --smoke
#   julia --project=. scripts/run_two_choice.jl 6 110000 20261101
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

# The entire stochastic environment, conditional forecasts, marginal utility
# continuation and realized household equations are the existing implementation.
# Only the firm's adjustment choice set differs.
function simulate_two_choice_period!(state::SimulationState)
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
        logmu = forecast.logmu_next + state.re.h * p.delta * epsilon
        mu = s.mu * exp(logmu)
        isfinite(mu) && mu > 0.0 ||
            throw(DomainError(mu, "invalid two-choice continuation marginal utility"))
        realized_macro = Nirei2015._realized_macroeconomy(
            state.calibration, kp, chosen_X, logmu, logmu_prior
        )
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
            capital_linear_gap=log(chosen_K / s.K) -
                               (forecast.logK_next + p.delta*epsilon),
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
    T > 0 && N > 0 && replication > 0 ||
        throw(ArgumentError("T, N and replication must be positive"))
    state = initialize_simulation(; N=N, base_seed=seed, replication=replication)
    records = NamedTuple[]
    failure = nothing
    for _ in 1:T
        try
            push!(records, simulate_two_choice_period!(state))
        catch err
            failure = sprint(showerror, err)
            break
        end
    end
    root = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(root)
    stem = "two_choice_rep" * lpad(string(replication), 2, '0')
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
            "not_published_2015_replication" => true,
            "firms" => N, "seed" => seed, "replication" => replication,
            "selection" => "original Section 3.4 ES1 over X>0",
            "continuation" => "current first-order Appendix A mu innovation",
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
        )
        r = last(records)
        summary["last_period"] = Dict(
            "K_next" => r.K_next, "X" => r.X, "n_fixed" => r.n_fixed,
        )
    end
    if length(records) > 101
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
    println("Two-choice counterfactual: T=$T N=$N seed=$seed replication=$replication")
    println("Completed ", length(records), "/", T, " periods")
    if !isempty(records)
        r = first(records)
        println("t=0: X=", r.X, " Knext=", r.K_next,
                " fixed points=", r.n_fixed, " positive-X=", r.n_positive,
                " investing=", r.n_investing)
        r = last(records)
        println("last: t=", r.t, " X=", r.X, " Knext=", r.K_next)
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
