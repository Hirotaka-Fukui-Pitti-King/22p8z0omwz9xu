#!/usr/bin/env julia
# Windows PowerShell: julia --project=. scripts/run_equilibrium.jl

using Nirei2015
using TOML
import Pkg

function main()
    calibration = run_calibration()
    re = compute_linear_re_solution(calibration)
    base_seed, replication = 20261101, 1
    rngs = firm_rngs(base_seed, replication)
    firms = initialize_firms(calibration, rngs)
    K_prior = aggregate_capital(firms, calibration)
    expectation = compute_saddle_path_expectations(re, calibration, K_prior)
    X_expected = expectation.X_next
    workspace = initialize_equilibrium_workspace(firms)

    draw_productivity!(firms, rngs, calibration)
    t0 = time_ns()
    solution = solve_equilibrium!(
        workspace, firms, calibration, X_expected, K_prior
    )
    elapsed_seconds = (time_ns() - t0) * 1.0e-9

    K_after = aggregate_capital(firms, calibration)
    abs(K_after - solution.K_next) <= 1e-10 * max(1.0, solution.K_next) ||
        error("Post-advance aggregate capital is inconsistent")
    !firms.next_ready || error("Next-productivity buffer not cleared after advancement")

    data = Dict(
        "source" => Dict(
            "paper" => "Nirei (2015), Theoretical Economics 10, 953-985",
            "doi" => "10.3982/TE1611",
            "equations" => "(10), (13), (14), (16); Section 3.4 ES1",
            "timing" => "K_prior is computed using a_i,t and k_i,t; equilibrium uses a_i,t+1",
            "selection" => "minimum absolute log investment innovation across positive-X fixed points",
            "tie_break" => "closest log K to continuum expected K, then smaller K",
            "nonpositive_X" => "excluded from log-based selection and counted (not clipped)",
            "intervals" => "(left,right], events at right applied after testing the interval",
            "expected_level" => "exp of first-order expected log deviation, not exact E[X]",
            "initialization" => "saddle-path mu for standalone one-period example only",
        ),
        "run" => Dict(
            "firms" => length(firms.capital),
            "base_seed" => base_seed,
            "replication" => replication,
            "rng_type" => "Random.Xoshiro",
            "elapsed_seconds" => elapsed_seconds,
        ),
        "equilibrium" => Dict(
            "K_prior" => solution.K_prior,
            "K_expected_continuum" => solution.expected_K,
            "X_expected" => solution.expected_X,
            "K_selected" => solution.K_next,
            "X_realized" => solution.X,
            "epsilon_log_investment" => solution.epsilon,
            "es1_objective" => solution.objective,
            "all_fixed_points" => solution.n_fixed,
            "admissible_positive_X" => solution.n_admissible,
            "excluded_nonpositive_X" => solution.n_nonpositive,
            "minimum_fixed_point_K" => solution.minimum_K,
            "maximum_fixed_point_K" => solution.maximum_K,
            "investing_firms" => solution.n_investing,
            "inactive_firms" => solution.n_inactive,
            "divesting_firms" => solution.n_divesting,
        ),
        "residuals" => Dict(
            "fixed_point_K" => solution.capital_residual,
            "aggregate_investment" => solution.investment_residual,
            "threshold_violation" => solution.max_threshold_violation,
            "post_advance_K" => K_after - solution.K_next,
        ),
        "environment" => Dict(
            "julia_version" => string(VERSION),
            "active_project" => something(Base.active_project(), ""),
            "dependencies" => Dict(
                info.name => (isnothing(info.version) ? "stdlib" : string(info.version))
                for info in values(Pkg.dependencies())
                if info.name in ("Distributions", "QuadGK", "TOML", "LinearAlgebra", "Random", "Statistics")
            ),
        ),
    )

    output = normpath(joinpath(@__DIR__, "..", "results", "equilibrium.toml"))
    mkpath(dirname(output))
    open(output, "w") do io
        TOML.print(io, data)
    end

    println("Nirei (2015), Issue #4 — finite-firm ES1")
    println("N = ", length(firms.capital), "; elapsed seconds = ", elapsed_seconds)
    println("Fixed points = ", solution.n_fixed,
            "; positive-X = ", solution.n_admissible,
            "; nonpositive-X = ", solution.n_nonpositive)
    println("Expected X = ", solution.expected_X, "; realized X = ", solution.X)
    println("Expected K* = ", solution.expected_K, "; selected K = ", solution.K_next)
    println("ES1 epsilon = ", solution.epsilon, "; objective = ", solution.objective)
    println("Fixed-point residual = ", solution.capital_residual,
            "; investment residual = ", solution.investment_residual)
    println("Result written: ", output)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
