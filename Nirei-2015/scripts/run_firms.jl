#!/usr/bin/env julia
# From Nirei-2015/: julia --project=. scripts/run_firms.jl

using Nirei2015
using Statistics: mean, std
using Distributions: Normal, cdf
using TOML
import Pkg

function sample_statistics(v::Vector{Float64})
    return Dict("mean" => mean(v), "std" => std(v),
                "minimum" => minimum(v), "maximum" => maximum(v))
end

function main()
    calibration = run_calibration()
    base_seed = 20261101
    replication = 1
    rngs = firm_rngs(base_seed, replication)
    firms = initialize_firms(calibration, rngs)
    n = length(firms.capital)
    p, c = calibration.parameters, calibration.constants

    # Independent fresh initialization with identical stream seeds.
    replay = initialize_firms(calibration, firm_rngs(base_seed, replication))
    deterministic_replay = firms.lambda == replay.lambda &&
                           firms.capital == replay.capital &&
                           firms.productivity == replay.productivity &&
                           firms.initial_gap == replay.initial_gap

    K_summary = realized_initial_gap_summary(firms, calibration)
    lambda_before = copy(firms.lambda)
    capital_before = copy(firms.capital)
    productivity_before = copy(firms.productivity)
    a_next = draw_productivity!(firms, rngs, calibration)
    unchanged_lambda = firms.lambda == lambda_before
    unchanged_capital = firms.capital == capital_before
    unchanged_productivity = firms.productivity == productivity_before
    next_draw_positive = all(x -> isfinite(x) && x > 0.0, a_next)
    initial_gap_valid = all(x -> 0.0 <= x < 1.0, firms.initial_gap)
    lambda_low = p.lambda_mean - p.lambda_truncation_sd * p.lambda_sd
    lambda_high = p.lambda_mean + p.lambda_truncation_sd * p.lambda_sd
    lambda_support_valid = all(x -> lambda_low <= x <= lambda_high, firms.lambda)
    all((deterministic_replay, unchanged_lambda, unchanged_capital,
         unchanged_productivity, next_draw_positive,
         initial_gap_valid, lambda_support_valid)) ||
        error("Firm-state validity / reproducibility check failed")

    log_initial_a = log.(firms.productivity)
    sample_a_power = mean(x -> x^c.productivity_exponent, firms.productivity)
    expected_lambda_mean = p.lambda_mean
    # Conditional variance of a symmetric ±d normal truncation.
    # Theoretical sample SD differs from p.lambda_sd due to truncation.
    z = p.lambda_truncation_sd
    norm_const = 2.0 * cdf(Normal(), z) - 1.0
    variance_factor = 1.0 - 2.0 * z * exp(-z^2 / 2.0) /
                                  (sqrt(2.0 * pi) * norm_const)
    theoretical_truncated_sd = p.lambda_sd * sqrt(variance_factor)
    seeds = stream_seeds(rngs)

    output = Dict(
        "source" => Dict(
            "paper" => "Nirei (2015), Sections 2.2, 2.6, 3.2-3.3, Appendix C",
            "initialization" => "Kbar-based sampled capital gaps; finite-N aggregate capital not rescaled",
            "note" => "Realized finite-N gaps relative to K0 may differ from originally drawn uniform gaps",
        ),
        "environment" => Dict(
            "julia_version" => string(VERSION),
            "active_project" => something(Base.active_project(), ""),
            "dependencies" => Dict(info.name => (isnothing(info.version) ? "stdlib" : string(info.version))
                for info in values(Pkg.dependencies())
                if info.name in ("Distributions", "QuadGK", "Random", "Statistics", "TOML")),
        ),
        "rng" => Dict(
            "rng_type" => "Random.Xoshiro",
            "base_seed" => base_seed,
            "replication" => replication,
            "seed_lambda" => string(seeds.lambda),
            "seed_initial_productivity" => string(seeds.initial_productivity),
            "seed_initial_gap" => string(seeds.initial_gap),
            "seed_subsequent_productivity" => string(seeds.subsequent_productivity),
        ),
        "parameters" => Dict(
            "firms" => n,
            "lambda_mean_untruncated" => p.lambda_mean,
            "lambda_sd_untruncated" => p.lambda_sd,
            "lambda_lower" => lambda_low,
            "lambda_upper" => lambda_high,
            "lambda_theoretical_mean_truncated" => expected_lambda_mean,
            "lambda_theoretical_sd_truncated" => theoretical_truncated_sd,
            "mean_log_productivity_theoretical" => c.mean_log_a,
            "sd_log_productivity_theoretical" => p.sigma_log_a,
            "E_a_power_theoretical" => 1.0,
            "productivity_exponent" => c.productivity_exponent,
        ),
        "sample_lambda" => sample_statistics(firms.lambda),
        "sample_log_productivity" => sample_statistics(log_initial_a),
        "sample_initial_gap" => sample_statistics(firms.initial_gap),
        "sample_initial_capital" => sample_statistics(firms.capital),
        "sample_next_productivity" => sample_statistics(a_next),
        "moments" => Dict(
            "sample_mean_a_power" => sample_a_power,
            "finite_K0" => K_summary.K0,
            "steady_Kbar" => K_summary.Kbar,
            "K0_relative_deviation" => K_summary.K0_relative_deviation,
            "realized_gap_mean" => K_summary.realized_gap_mean,
            "realized_gap_std" => K_summary.realized_gap_std,
            "realized_gap_min" => K_summary.realized_gap_min,
            "realized_gap_max" => K_summary.realized_gap_max,
            "realized_gap_outside_0_1_count" => K_summary.realized_gap_outside_0_1_count,
        ),
        "checks" => Dict(
            "same_seed_initialization_bitwise" => deterministic_replay,
            "lambda_unchanged_by_productivity_draw" => unchanged_lambda,
            "capital_unchanged_by_productivity_draw" => unchanged_capital,
            "current_productivity_unchanged_by_next_draw" => unchanged_productivity,
            "next_productivity_positive_finite" => next_draw_positive,
            "sampled_gaps_in_0_1" => initial_gap_valid,
            "lambda_in_truncated_support" => lambda_support_valid,
        ),
    )

    path = normpath(joinpath(@__DIR__, "..", "results", "firms.toml"))
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, output)
    end
    println("Nirei (2015): Issue #3 firm heterogeneity and shocks")
    println("Julia version: ", VERSION)
    println("N = ", n, ", RNG = Xoshiro, seed = ", base_seed)
    println("Sample lambda mean = ", mean(firms.lambda),
            "; SD = ", std(firms.lambda),
            "; theoretical truncated SD = ", theoretical_truncated_sd)
    println("Sample E[a^q] = ", sample_a_power)
    println("K0 = ", K_summary.K0, "; Kbar = ", K_summary.Kbar,
            "; K0/Kbar - 1 = ", K_summary.K0_relative_deviation)
    println("Finite-N realized gaps outside [0,1): ",
            K_summary.realized_gap_outside_0_1_count)
    println("Deterministic replay: ", deterministic_replay)
    println("Result written: ", path)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
