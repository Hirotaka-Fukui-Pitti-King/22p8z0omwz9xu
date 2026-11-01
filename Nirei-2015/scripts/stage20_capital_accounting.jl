#!/usr/bin/env julia
# Stage 20: exact CES growth accounting and stationary-gap audit.
# The original firm policy, ES1, macro continuation and RNG are UNCHANGED.
# Commands from Nirei-2015:
# julia --project=. scripts/stage20_capital_accounting.jl --smoke
# julia --project=. scripts/stage20_capital_accounting.jl 20 110000 20261101
# julia --project=. scripts/stage20_capital_accounting.jl 1100 110000 20261101
using Nirei2015
using TOML
using Statistics: mean
include(joinpath(@__DIR__, "run_two_choice.jl"))

# This verifies an identity, not an approximation or a chosen economic closure.
function growth_identity(aold, anew, kold, knew, alpha, rho, delta)
    n = length(kold)
    length(aold) == n && length(anew) == n && length(knew) == n ||
        throw(DimensionMismatch("firm-state arrays must agree"))
    S_old = 0.0
    S_new_no_adjust = 0.0
    S_new_adjust = 0.0
    gross = inv(1.0 - delta)
    @inbounds for i in eachindex(kold)
        old_power = kold[i]^rho
        old_weight = aold[i]^(rho / alpha)
        new_weight = anew[i]^(rho / alpha)
        S_old += old_weight * old_power
        S_new_no_adjust += new_weight * old_power
        S_new_adjust += new_weight * (knew[i] * gross)^rho
    end
    depreciation = log1p(-delta)
    productivity = log(S_new_no_adjust / S_old) / rho
    adjustment = log(S_new_adjust / S_new_no_adjust) / rho
    K_old = (S_old / n)^(inv(rho))
    K_new = (S_new_adjust / n)^(inv(rho)) * (1.0 - delta)
    actual = log(K_new / K_old)
    residual = actual - depreciation - productivity - adjustment
    abs(residual) <= 1e-11 ||
        throw(ErrorException("exact CES growth decomposition failed: $residual"))
    return (
        depreciation=depreciation, productivity=productivity,
        adjustment=adjustment, growth=actual, residual=residual,
        K_prior_check=K_old, K_next_check=K_new,
    )
end

# Reconstruct the exact gap s=log(k/(b(a)*K))/log(lambda).
# Check both support of s and closeness of fractional gaps to uniform.
# The grouped CDF measure is a diagnostic, not a statistical significance test.
function gap_summary(firms::FirmState, state::SimulationState, oldcapital, olda, K_before)
    n = length(firms.capital)
    buckets = zeros(Int, 10)
    outside = 0
    total = 0.0
    total_sq = 0.0
    low = Inf
    high = -Inf
    @inbounds for i in 1:n
        lambda = firms.lambda[i]
        s = log(oldcapital[i] /
                (threshold_coefficient(lambda, olda[i], state.calibration) * K_before)) /
            log(lambda)
        total += s
        total_sq += s*s
        low = min(low, s)
        high = max(high, s)
        if s < 0.0 || s >= 1.0
            outside += 1
        end
        f = s - floor(s)
        index = min(10, max(1, floor(Int, 10*f) + 1))
        buckets[index] += 1
    end
    cumsum_count = 0
    cdf_error = 0.0
    for j in 1:10
        cumsum_count += buckets[j]
        cdf_error = max(cdf_error, abs(cumsum_count/n - j/10))
    end
    avg = total / n
    return (
        gap_mean=avg, gap_std=sqrt(max(0.0, total_sq/n - avg^2)),
        gap_min=low, gap_max=high, outside_share=outside/n,
        frac_gap_cdf_decile_error=cdf_error,
    )
end

# Reaction to a fixed candidate aggregate capital, holding the same
# next-productivity shocks and existing pre-investment firm stocks.
function conditional_reaction(firms, calibration, oldcapital, anew, candidateK)
    p = calibration.parameters
    rho = calibration.constants.rho
    n = length(oldcapital)
    power = 0.0
    investment = 0.0
    n_invest = 0
    @inbounds for i in 1:n
        d = (1.0 - p.delta) * oldcapital[i]
        b = threshold_coefficient(firms.lambda[i], anew[i], calibration)
        invest = candidateK > d / b
        ki = invest ? d*firms.lambda[i] : d
        power += anew[i]^(rho / p.alpha) * ki^rho
        investment += ki - d
        n_invest += invest
    end
    return (
        responseK=(power/n)^(inv(rho)),
        responseX=investment/n,
        investing_share=n_invest/n,
    )
end

function collect_accounting(state, record)
    # After simulate_two_choice_period!, the O(1) swaps retain the exact
    # PREVIOUS stocks in workspace.next_capital and next_productivity.
    f = state.firms
    oldk = state.workspace.next_capital
    olda = f.next_productivity
    newk = f.capital
    anew = f.productivity
    p = state.calibration.parameters
    rho = state.calibration.constants.rho
    terms = growth_identity(olda, anew, oldk, newk, p.alpha, rho, p.delta)
    abs(terms.K_prior_check - record.K_prior) < 1e-9 * max(1.0,record.K_prior) ||
        throw(ErrorException("old firm arrays no longer match K_prior"))
    abs(terms.K_next_check - record.K_next) < 1e-9 * max(1.0,record.K_next) ||
        throw(ErrorException("new firm arrays no longer match K_next"))
    oldg = gap_summary(f, state, oldk, olda, record.K_prior)
    nextg = gap_summary(f, state, newk, anew, record.K_next)
    expectation_reaction = conditional_reaction(
        f, state.calibration, oldk, anew, record.expected_K
    )
    prior_reaction = conditional_reaction(
        f, state.calibration, oldk, anew, record.K_prior
    )
    return (
        t=record.t, K_prior=record.K_prior, K_next=record.K_next,
        X=record.X, expected_X=record.expected_X,
        selected_investment_share=record.n_investing/length(newk),
        investment_ratio=record.X/record.expected_X,
        fixed_points=record.n_fixed, stationary_X=state.calibration.steady_state.X,
        depreciation=terms.depreciation, productivity=terms.productivity,
        adjustment=terms.adjustment, growth=terms.growth,
        identity_residual=terms.residual,
        old_gap_mean=oldg.gap_mean, old_gap_std=oldg.gap_std,
        old_gap_outside_share=oldg.outside_share,
        new_gap_mean=nextg.gap_mean, new_gap_std=nextg.gap_std,
        new_gap_min=nextg.gap_min, new_gap_max=nextg.gap_max,
        new_gap_outside_share=nextg.outside_share,
        new_fractional_gap_cdf_decile_error=nextg.frac_gap_cdf_decile_error,
        expectedK_reactionK=expectation_reaction.responseK,
        expectedK_reactionX=expectation_reaction.responseX,
        expectedK_reaction_residual=expectation_reaction.responseK-record.expected_K,
        priorK_reactionK=prior_reaction.responseK,
        priorK_reaction_residual=prior_reaction.responseK-record.K_prior,
    )
end

function accounting_smoke()
    alpha, rho, delta = 0.26, 0.49, 0.02
    olda = [0.98, 1.02, 1.01]
    anew = [1.01, 0.99, 1.00]
    oldk = [1.2, 1.8, 0.7]
    lambda = [1.028, 1.021, 1.036]
    choice = [1, 0, 1]
    newk = [(1.0-delta)*oldk[i]*lambda[i]^choice[i] for i in 1:3]
    result = growth_identity(olda, anew, oldk, newk, alpha, rho, delta)
    abs(result.residual) < 1e-12 || error("independent CES identity smoke failed")
    # Root set scales exactly with all old firm capital for any positive q.
    state = initialize_simulation(; N=8, base_seed=20261101)
    draw_productivity!(state.firms, state.rngs, state.calibration)
    roots = enumerate_two_choice(state.firms, state.calibration)
    q = 1.37
    f = state.firms
    scaled = FirmState(copy(f.lambda), q .* f.capital,
                       copy(f.productivity), copy(f.next_productivity),
                       copy(f.initial_gap), copy(f.initial_b), true)
    scaled_roots = enumerate_two_choice(scaled, state.calibration)
    length(roots.K) == length(scaled_roots.K) ||
        error("capital homogeneity root count mismatch")
    for i in eachindex(roots.K)
        isapprox(scaled_roots.K[i], q*roots.K[i]; rtol=1e-11) ||
            error("homogeneity failed for capital root")
        isapprox(scaled_roots.X[i], q*roots.X[i]; rtol=1e-11, atol=1e-12) ||
            error("homogeneity failed for investment root")
    end
    println("PASS: independent exact CES growth identity; homogeneous root-set scaling")
end

function write_table(path, records)
    isempty(records) && return
    keys = propertynames(first(records))
    open(path, "w") do io
        println(io, join(string.(keys), ","))
        for record in records
            println(io, join((string(getproperty(record, k)) for k in keys), ","))
        end
    end
end

function main_stage20()
    if !isempty(ARGS) && ARGS[1] == "--smoke"
        accounting_smoke()
        return nothing
    end
    T = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
    N = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 110000
    seed = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 20261101
    T > 0 && N > 0 || throw(ArgumentError("T,N must be positive"))
    state = initialize_simulation(; N=N, base_seed=seed, replication=1)
    records = NamedTuple[]
    failure = nothing
    for _ in 1:T
        try
            chosen = simulate_two_choice_period!(state, :linear)
            push!(records, collect_accounting(state, chosen))
        catch err
            failure = sprint(showerror, err)
            break
        end
    end
    output_root = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(output_root)
    stem = joinpath(output_root, "stage20_capital_accounting_T$(T)")
    write_table(stem*".csv", records)
    stat = Dict{String,Any}(
        "completed" => length(records) == T && failure === nothing,
        "completed_periods" => length(records),
        "requested_periods" => T,
        "failure" => something(failure, "none"),
        "N" => N, "seed" => seed,
    )
    if !isempty(records)
        stat["maximum_growth_identity_error"] =
            maximum(r -> abs(r.identity_residual), records)
        stat["whole_path_mean_components"] = Dict(
            "depreciation" => mean(r.depreciation for r in records),
            "productivity" => mean(r.productivity for r in records),
            "adjustment" => mean(r.adjustment for r in records),
            "capital_growth" => mean(r.growth for r in records),
        )
        for (label, selected) in (
            ("initial_100", [r for r in records if r.t < 100]),
            ("after_100", [r for r in records if r.t >= 100]),
        )
            if !isempty(selected)
                stat[label] = Dict(
                    "periods" => length(selected),
                    "mean_depreciation" => mean(r.depreciation for r in selected),
                    "mean_productivity" => mean(r.productivity for r in selected),
                    "mean_adjustment" => mean(r.adjustment for r in selected),
                    "mean_growth" => mean(r.growth for r in selected),
                    "mean_gap_outside_share" => mean(r.new_gap_outside_share for r in selected),
                    "mean_gap_cdf_error" =>
                        mean(r.new_fractional_gap_cdf_decile_error for r in selected),
                )
            end
        end
    end
    open(stem*"_summary.toml", "w") do io
        TOML.print(io, stat)
    end
    println("Stage 20 exact CES capital accounting: completed ", length(records), "/", T)
    if !isempty(records)
        for t in (0,5,19,139,1099)
            if t < length(records)
                r = records[t+1]
                println("t=",t," K=",r.K_next," X=",r.X,
                        " log_growth=",r.growth,
                        " depreciation=",r.depreciation,
                        " productivity=",r.productivity,
                        " adjustment=",r.adjustment,
                        " expectedK_reaction_residual=",r.expectedK_reaction_residual,
                        " gap_outside_share=",r.new_gap_outside_share)
            end
        end
        println("Mean components: ", stat["whole_path_mean_components"])
    end
    println("CSV: ", stem*".csv")
    println("Summary: ", stem*"_summary.toml")
    failure === nothing || error("Stage 20 failed: "*failure)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_stage20()
end
