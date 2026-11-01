#!/usr/bin/env julia
# Original 2015 three-choice benchmark, UNCHANGED simulate_period! / ES1.
# Screen an ex-ante fixed consecutive set of independently seeded paths.
# Windows PowerShell from Nirei-2015:
#   julia --project=. scripts/scan_original_seeds.jl 6 110000 20261102 20
# Arguments: periods, firms, first_seed, number_of_seeds.
# Failures are recorded (not filtered, resampled, or silently retried).
using Nirei2015
using TOML
using Statistics: std

function csv_field(x)
    x === nothing && return ""
    raw = replace(string(x), "\r" => " ", "\n" => " ")
    return "\"" * replace(raw, "\"" => "\"\"") * "\""
end

function collect_one_seed(calibration::CalibrationResult, T::Int, N::Int, seed::Int)
    state = initialize_simulation(; calibration=calibration, N=N,
                                   base_seed=seed, replication=1)
    initial_K = state.K
    records = PeriodRecord[]
    sizehint!(records, T)
    failure = "none"
    failure_class = "none"
    start_time = time_ns()
    for _ in 1:T
        try
            push!(records, simulate_period!(state))
        catch err
            err isa InterruptException && rethrow()
            failure = sprint(showerror, err)
            failure_class = occursin("No equilibrium has X>0", failure) ?
                "no_positive_X_equilibrium" : string(typeof(err))
            break
        end
    end
    elapsed = (time_ns() - start_time) * 1.0e-9
    n_completed = length(records)
    completed = n_completed == T
    max_capital_residual = isempty(records) ? nothing :
        maximum(r -> abs(r.fixed_point_residual), records)
    max_investment_residual = isempty(records) ? nothing :
        maximum(r -> abs(r.investment_residual), records)
    max_level_resource_gap = isempty(records) ? nothing :
        maximum(r -> abs(r.level_resources_gap), records)
    max_abs_epsilon = isempty(records) ? nothing :
        maximum(r -> abs(r.epsilon), records)
    # As in the benchmark, burn-in = 100; report moments only if >= 2 retained observations.
    retained = n_completed > 101 ? @view(records[101:end]) : PeriodRecord[]
    percent_std(field) = length(retained) >= 2 ?
        100 * std(getproperty(r, field) for r in retained) : nothing
    return (
        seed=seed, requested_periods=T, completed_periods=n_completed,
        completed=completed, failure_period=completed ? nothing : state.period,
        failure_class=failure_class, initial_K=initial_K, final_K=state.K,
        first_X=isempty(records) ? nothing : first(records).X,
        last_X=isempty(records) ? nothing : last(records).X,
        max_capital_residual=max_capital_residual,
        max_investment_residual=max_investment_residual,
        max_level_resource_gap=max_level_resource_gap,
        max_abs_epsilon=max_abs_epsilon,
        Y_sd_percent=percent_std(:logY),
        C_sd_percent=percent_std(:logC),
        X_sd_percent=percent_std(:logX),
        L_sd_percent=percent_std(:logL),
        K_sd_percent=percent_std(:logK),
        elapsed_seconds=elapsed,
        failure=failure,
    )
end

function main()
    length(ARGS) == 4 || error(
        "Usage: julia --project=. scripts/scan_original_seeds.jl PERIODS FIRMS FIRST_SEED COUNT")
    T = parse(Int, ARGS[1])
    N = parse(Int, ARGS[2])
    first_seed = parse(Int, ARGS[3])
    count = parse(Int, ARGS[4])
    T >= 2 || throw(ArgumentError("PERIODS must be >= 2"))
    N >= 1 || throw(ArgumentError("FIRMS must be positive"))
    first_seed >= 1 || throw(ArgumentError("FIRST_SEED must be positive"))
    count >= 1 || throw(ArgumentError("COUNT must be positive"))
    count <= typemax(Int) - first_seed + 1 ||
        throw(ArgumentError("seed range exceeds Int"))
    cal = run_calibration()
    cal.parameters.sigma == 1.5 ||
        error("Expected original Table 1 risk aversion sigma=1.5")
    cal.parameters.firms == 110000 ||
        error("Unexpected original published number of firms")
    output_dir = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(output_dir)
    filename = "original_seed_scan_T$(T)_N$(N)_start$(first_seed)_count$(count)"
    csv_path = joinpath(output_dir, filename * ".csv")
    toml_path = joinpath(output_dir, filename * ".toml")
    n_success = 0
    n_no_positive = 0
    best_prefix = 0
    success_seeds = Int[]
    failures = Dict{String, Int}()
    fields = (
        :seed, :requested_periods, :completed_periods, :completed,
        :failure_period, :failure_class, :initial_K, :final_K,
        :first_X, :last_X, :max_capital_residual, :max_investment_residual,
        :max_level_resource_gap, :max_abs_epsilon,
        :Y_sd_percent, :C_sd_percent, :X_sd_percent, :L_sd_percent,
        :K_sd_percent, :elapsed_seconds, :failure,
    )
    println("Original published 2015 three-choice seed sensitivity; sigma=1.5, ES1 unchanged.")
    println("No shock rejection, no investment floor, no change to simulate_period!.")
    println("Seeds: $(first_seed) through $(first_seed + count - 1); T=$T, N=$N")
    # Stream CSV to disk after every seed; interrupted screening retains prior completed rows.
    open(csv_path, "w") do io
        println(io, join(string.(fields), ","))
        flush(io)
        for offset in 0:(count - 1)
            seed = first_seed + offset
            row = collect_one_seed(cal, T, N, seed)
            println(io, join((csv_field(getproperty(row, k)) for k in fields), ","))
            flush(io)
            best_prefix = max(best_prefix, row.completed_periods)
            if row.completed
                n_success += 1
                push!(success_seeds, seed)
            else
                failures[row.failure_class] = get(failures, row.failure_class, 0) + 1
                n_no_positive += row.failure_class == "no_positive_X_equilibrium"
            end
            println("seed=$seed periods=$(row.completed_periods)/$T " *
                    "status=$(row.completed ? "PASS" : row.failure_class) " *
                    "K_last=$(row.final_K) elapsed=$(round(row.elapsed_seconds; digits=2))s")
            flush(stdout)
        end
    end
    summary = Dict(
        "source" => "Original 2015 three-choice finite-N model (published sigma=1.5)",
        "method" => "Unmodified simulate_period!, original ES1, iid shocks; no resampling",
        "first_seed" => first_seed, "count" => count,
        "requested_periods" => T, "firms" => N,
        "n_completed" => n_success,
        "fraction_completed" => n_success / count,
        "n_no_positive_X" => n_no_positive,
        "longest_completed_prefix" => best_prefix,
        "success_seeds" => success_seeds,
        "failure_counts" => failures,
        "warning" => "Seed-conditional path completion does not establish global equilibrium, unconditional Table 2 replication, or nonlinear rational expectations."
    )
    open(toml_path, "w") do io
        TOML.print(io, summary)
    end
    println("Completed $n_success/$count paths through T=$T; " *
            "no-positive-X failures=$n_no_positive; longest prefix=$best_prefix")
    println("CSV: $csv_path")
    println("TOML: $toml_path")
    println("Any PASS seed is a CANDIDATE for T20, not a published benchmark replication.")
end

main()
