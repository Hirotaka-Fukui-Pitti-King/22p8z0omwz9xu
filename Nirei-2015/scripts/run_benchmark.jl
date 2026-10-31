#!/usr/bin/env julia
# From Nirei-2015/: julia --project=. scripts/run_benchmark.jl [T] [N] [seed]
# Defaults to T=1100, N=110000, seed=20261101. One trial only, Issue #5.

using Nirei2015
using TOML
using Statistics: mean, std
import Pkg

function csv_header()
    return string.(fieldnames(PeriodRecord))
end

function write_records_csv(path::String, records::Vector{PeriodRecord})
    open(path, "w") do io
        println(io, join(csv_header(), ","))
        for r in records
            println(io, join((string(getfield(r, key)) for key in fieldnames(PeriodRecord)), ","))
        end
    end
end

function run()
    T = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1100
    N = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 110000
    seed = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 20261101
    T >= 2 || throw(ArgumentError("T must be >= 2"))
    N >= 1 || throw(ArgumentError("N must be positive"))
    burnin = min(100, T - 1)
    state = initialize_simulation(; N=N, base_seed=seed, replication=1)
    initial_K = state.K
    initial_mu = state.mu_prior
    root = normpath(joinpath(@__DIR__, "..", "results"))
    mkpath(root)
    csv_path = joinpath(root, "simulation_rep01.csv")
    summary_path = joinpath(root, "simulation_summary.toml")

    start = time_ns()
    records = PeriodRecord[]
    fail = nothing
    try
        # Preserve the successful prefix for troubleshooting; no fake fill-in.
        sizehint!(records, T)
        for _ in 1:T
            push!(records, simulate_period!(state))
        end
    catch err
        fail = sprint(showerror, err)
    end
    elapsed = (time_ns() - start) * 1e-9
    write_records_csv(csv_path, records)

    summary = Dict(
        "status" => Dict(
            "completed" => fail === nothing && length(records) == T,
            "completed_periods" => length(records),
            "requested_periods" => T,
            "failure" => something(fail, "none"),
            "elapsed_seconds" => elapsed,
        ),
        "configuration" => Dict(
            "firms" => N, "periods" => T, "burnin" => burnin,
            "base_seed" => seed, "replication" => 1, "rng" => "Random.Xoshiro",
            "initial_K" => initial_K, "initial_mu" => initial_mu,
            "moment_calculation" => "Not part of Issue #5. Diagnostic summary only.",
            "conditional_expectations" => "Appendix A B^-1D two-state forecast of log deviations",
            "macro_closure" => "First-order bounded-continuation mu_t = h*log(K_{t+1}^firm/Kbar); finite-N forecast drift is measured",
            "price_setting" => "Predetermined price equations imposed in expectations, not realized macro identities",
            "capital" => "Firm aggregation is exact; Appendix A linear capital is a separate diagnostic",
        ),
        "environment" => Dict(
            "julia_version" => string(VERSION),
            "active_project" => something(Base.active_project(), ""),
            "packages" => Dict(info.name => (isnothing(info.version) ? "stdlib" : string(info.version))
                for info in values(Pkg.dependencies())
                if info.name in ("QuadGK", "Distributions", "Statistics", "Random", "TOML")),
        ),
        "integrity" => Dict(
            "max_abs_fixed_point_residual" => isempty(records) ? 0.0 : maximum(r -> abs(r.fixed_point_residual), records),
            "max_abs_investment_residual" => isempty(records) ? 0.0 : maximum(r -> abs(r.investment_residual), records),
            "max_abs_production_linear_residual" => isempty(records) ? 0.0 : maximum(r -> abs(r.production_linear_residual), records),
            "max_abs_resources_linear_residual" => isempty(records) ? 0.0 : maximum(r -> abs(r.resources_linear_residual), records),
            "max_abs_marginal_utility_linear_residual" => isempty(records) ? 0.0 : maximum(r -> abs(r.marginal_utility_linear_residual), records),
            "max_abs_stable_jump_residual" => isempty(records) ? 0.0 : maximum(r -> abs(r.stable_jump_residual), records),
        ),
    )
    if !isempty(records)
        summary["initial_period"] = Dict(
            "epsilon" => first(records).epsilon,
            "expected_X" => first(records).expected_X,
            "realized_X" => first(records).X,
            "K_next" => first(records).K_next_firm,
        )
    end
    if length(records) > burnin
        retained = @view records[(burnin+1):end]
        xs = [r.epsilon for r in retained]
        kls = [r.logK for r in retained]
        summary["post_burnin"] = Dict(
            "retained_periods" => length(retained),
            "mean_epsilon" => mean(xs),
            "sd_epsilon" => length(xs) > 1 ? std(xs) : 0.0,
            "mean_abs_epsilon" => mean(abs, xs),
            "maximum_abs_epsilon" => maximum(abs, xs),
            "mean_X" => mean(r.X for r in retained),
            "mean_expected_X" => mean(r.expected_X for r in retained),
            "mean_ES1_fixed_points" => mean(r.n_fixed for r in retained),
            "mean_investing_share" => mean(r.n_investing for r in retained)/N,
            "mean_capital_linear_gap" => mean(r.capital_linear_gap for r in retained),
            "final_K" => last(retained).K_next_firm,
            "logK_endpoint_trend_per_period" => length(retained)>1 ?
                (last(kls)-first(kls))/(length(kls)-1) : 0.0,
        )
    end

    open(summary_path, "w") do io
        TOML.print(io, summary)
    end
    println("Nirei (2015) Issue #5 dynamic simulation")
    println("T = ", T, "; N = ", N, "; seed = ", seed,
            "; completed periods = ", length(records), "; elapsed seconds = ", elapsed)
    println("CSV written: ", csv_path)
    println("Summary written: ", summary_path)
    if fail !== nothing
        println(stderr, "FAILED: ", fail)
        error("Simulation stopped at t=$(state.period); see summary and partial CSV")
    end
    println("Mean post-burnin epsilon: ", summary["post_burnin"]["mean_epsilon"])
    println("Final aggregate K: ", state.K)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
