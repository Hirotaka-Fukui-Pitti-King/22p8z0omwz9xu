#!/usr/bin/env julia
# Source-grounded historical alternative: Nirei 2013 binary one-sided investment.
# NOT the published 2015 three-choice policy. No src/ modifications.
# Run from Nirei-2015: julia --project=. scripts/historical_binary_protocol.jl [T] [N] [seed]
using Nirei2015
using TOML
using Statistics: std, cor

function binary_equilibrium!(state::SimulationState, expected_X::Float64)
    firms = state.firms
    cal = state.calibration
    p = cal.parameters
    rho = cal.constants.rho
    n = length(firms.capital)
    invn = inv(Float64(n))
    firms.next_ready || error("next productivity not drawn")
    switching = Vector{Float64}(undef, n)
    base_power = Vector{Float64}(undef, n)
    increments = Vector{Float64}(undef, n)
    depreciated = Vector{Float64}(undef, n)
    base_sum = 0.0
    for i in 1:n
        d = (1.0 - p.delta) * firms.capital[i]
        a = firms.next_productivity[i]
        lambda = firms.lambda[i]
        b = threshold_coefficient(lambda, a, cal)
        switching[i] = d / b
        base_power[i] = a^(rho / p.alpha) * d^rho * invn
        increments[i] = (lambda - 1.0) * d * invn
        depreciated[i] = d
        base_sum += base_power[i]
    end
    order = sortperm(switching)
    candidates_K = Float64[]
    candidates_X = Float64[]
    power = base_sum
    investment = 0.0
    previous = 0.0
    j = 1
    while j <= n
        boundary = switching[order[j]]
        K = power^(inv(rho))
        if previous < K <= boundary
            push!(candidates_K, K)
            push!(candidates_X, investment)
        end
        while j <= n && switching[order[j]] == boundary
            i = order[j]
            power += base_power[i] * (firms.lambda[i]^rho - 1.0)
            investment += increments[i]
            j += 1
        end
        previous = boundary
    end
    K = power^(inv(rho))
    if K > previous
        push!(candidates_K, K)
        push!(candidates_X, investment)
    end
    best = 0
    criterion = Inf
    for i in eachindex(candidates_K)
        x = candidates_X[i]
        if x > 0.0
            dist = abs(log(x / expected_X))
            if dist < criterion
                best = i
                criterion = dist
            end
        end
    end
    best > 0 || error("historical binary rule has no positive-X equilibrium; roots=$(length(candidates_K))")
    selected_K = candidates_K[best]
    selected_X = candidates_X[best]
    count_invest = 0
    next_k = Vector{Float64}(undef, n)
    sum_x = 0.0
    for i in 1:n
        invest = selected_K > switching[i]
        next_k[i] = invest ? depreciated[i] * firms.lambda[i] : depreciated[i]
        count_invest += invest
        sum_x += next_k[i] - depreciated[i]
    end
    exact_K = aggregate_capital(firms.next_productivity, next_k, p.alpha, rho)
    exact_X = sum_x * invn
    residual_K = abs(exact_K - selected_K)
    residual_X = abs(exact_X - selected_X)
    residual_K <= 1e-9 * max(1.0, abs(selected_K)) ||
        error("direct CES reaggregation failed: $residual_K")
    residual_X <= 1e-9 * max(1.0, abs(selected_X)) ||
        error("direct investment reaggregation failed: $residual_X")
    firms.capital = next_k
    advance_productivity!(firms)
    return (K=selected_K, X=selected_X, count_invest=count_invest,
            roots=length(candidates_K), criterion=criterion,
            residual_K=residual_K, residual_X=residual_X)
end

function binary_period!(state::SimulationState)
    ss = state.calibration.steady_state
    t = state.period
    kp = state.K
    mup = state.mu_prior
    ex = compute_expected_equilibrium(state.re, state.calibration, kp, mup)
    draw_productivity!(state.firms, state.rngs, state.calibration)
    eq = binary_equilibrium!(state, ex.X_next)
    logK_prior = log(kp / ss.K)
    logmu_prior = log(mup / ss.mu)
    logK_next = log(eq.K / ss.K)
    logmu = state.re.h * logK_next
    mu = ss.mu * exp(logmu)
    logX = log(eq.X / ss.X)
    vals = Nirei2015._realized_macroeconomy(
        state.calibration, state.re, logK_prior, logX, logmu, logmu_prior)
    resource_gap = (vals.Y - vals.C - eq.X) / vals.Y
    isfinite(resource_gap) || error("nonfinite level resource gap")
    state.K = eq.K
    state.mu_prior = mu
    state.period += 1
    return (t=t, K=eq.K, X=eq.X, expected_X=ex.X_next,
            epsilon=log(eq.X / ex.X_next), roots=eq.roots,
            investing=eq.count_invest, Y=vals.Y, C=vals.C, L=vals.L,
            w=vals.w, mu=mu, residual_K=eq.residual_K,
            residual_X=eq.residual_X, resource_gap=resource_gap)
end

function main()
    T = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
    N = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 110000
    seed = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 20261101
    T >= 1 && N >= 1 || error("T and N must be positive")
    state = initialize_simulation(N=N, base_seed=seed)
    out = joinpath(@__DIR__, "..", "results")
    mkpath(out)
    records = NamedTuple[]
    failure = "none"
    for _ in 1:T
        try
            r = binary_period!(state)
            push!(records, r)
            if r.t <= 5 || r.t == T - 1
                println("t=$(r.t) K=$(r.K) X=$(r.X) Xe=$(r.expected_X) epsilon=$(r.epsilon) roots=$(r.roots) investing=$(r.investing)")
            end
        catch e
            failure = sprint(showerror, e)
            println(stderr, "HISTORICAL BINARY FAILURE at t=$(state.period): $failure")
            break
        end
    end
    prefix = "historical_binary_T$(T)_N$(N)_seed$(seed)"
    csvpath = joinpath(out, prefix * ".csv")
    open(csvpath, "w") do io
        headers = (:t,:K,:X,:expected_X,:epsilon,:roots,:investing,:Y,:C,:L,:w,:mu,
                   :residual_K,:residual_X,:resource_gap)
        println(io, join(string.(headers), ","))
        for r in records
            println(io, join((string(getfield(r, h)) for h in headers), ","))
        end
    end
    completed = length(records) == T
    summary = Dict(
        "variant" => "Nirei 2013 binary one-sided choices, 2015 calibration; NOT 2015 three-choice benchmark",
        "T" => T, "N" => N, "seed" => seed, "completed" => completed,
        "completed_periods" => length(records), "failure" => failure,
        "max_abs_CES_residual" => isempty(records) ? 0.0 : maximum(r -> r.residual_K, records),
        "max_abs_investment_residual" => isempty(records) ? 0.0 : maximum(r -> r.residual_X, records),
        "max_abs_level_resource_gap" => isempty(records) ? 0.0 : maximum(r -> abs(r.resource_gap), records),
    )
    if length(records) > 100
        kept = records[101:end]
        summary["postburnin"] = Dict(
            "Y_sd_percent" => 100 * std(log(r.Y) for r in kept),
            "C_sd_percent" => 100 * std(log(r.C) for r in kept),
            "X_sd_percent" => 100 * std(log(r.X) for r in kept),
            "L_sd_percent" => 100 * std(log(r.L) for r in kept),
            "K_sd_percent" => 100 * std(log(r.K) for r in kept),
            "Y_C_correlation" => cor([log(r.Y) for r in kept], [log(r.C) for r in kept]),
            "Y_X_correlation" => cor([log(r.Y) for r in kept], [log(r.X) for r in kept]),
            "Y_L_correlation" => cor([log(r.Y) for r in kept], [log(r.L) for r in kept]),
            "Y_K_correlation" => cor([log(r.Y) for r in kept], [log(r.K) for r in kept]),
        )
    end
    tomlpath = joinpath(out, prefix * ".toml")
    open(tomlpath, "w") do io
        TOML.print(io, summary)
    end
    println("RESULT: ", summary)
    println("Output: ", csvpath, " / ", tomlpath)
    completed || error("historical binary test incomplete at t=$(state.period)")
end

main()
