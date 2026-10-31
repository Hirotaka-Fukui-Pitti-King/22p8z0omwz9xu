#!/usr/bin/env julia
# Run from Nirei-2015/: julia --project=. scripts/run_calibration.jl

using Nirei2015
using TOML
import Pkg

function parameter_dictionary(p::BenchmarkParameters)
    return Dict(
        "alpha" => p.alpha,
        "delta" => p.delta,
        "beta" => p.beta,
        "sigma" => p.sigma,
        "eta" => p.eta,
        "psi" => p.psi,
        "zeta" => p.zeta,
        "lambda_mean" => p.lambda_mean,
        "lambda_sd_before_truncation" => p.lambda_sd,
        "lambda_truncation_sd" => p.lambda_truncation_sd,
        "lambda_lower" => p.lambda_mean - p.lambda_truncation_sd * p.lambda_sd,
        "lambda_upper" => p.lambda_mean + p.lambda_truncation_sd * p.lambda_sd,
        "sigma_log_a" => p.sigma_log_a,
        "firms" => p.firms,
    )
end

function dependency_versions()
    versions = Dict{String,String}()
    for info in values(Pkg.dependencies())
        if info.name in ("Distributions", "QuadGK", "TOML")
            versions[info.name] = isnothing(info.version) ? "stdlib" : string(info.version)
        end
    end
    return versions
end

function save_calibration(path::AbstractString, result::CalibrationResult)
    p, c, a, s = result.parameters, result.constants, result.threshold, result.steady_state
    residuals = equilibrium_residuals(result)

    data = Dict(
        "source" => Dict(
            "paper" => "Nirei (2015), An Interaction-Based Foundation of Aggregate Investment Fluctuations",
            "doi" => "10.3982/TE1611",
            "equations" => "Table 1, Appendix C, (7)-(12), (16)-(18)",
            "status" => "computed locally; not a published table of A or steady-state levels",
        ),
        "environment" => Dict(
            "julia_version" => string(VERSION),
            "dependencies" => dependency_versions(),
            "active_project" => something(Base.active_project(), ""),
        ),
        "parameters" => parameter_dictionary(p),
        "structural_constants" => Dict(
            "c_L" => c.c_L,
            "rho" => c.rho,
            "productivity_exponent" => c.productivity_exponent,
            "mean_log_a" => c.mean_log_a,
        ),
        "threshold_constant" => Dict(
            "A" => a.A,
            "integral" => a.integral,
            "integral_error_estimate" => a.integral_error,
        ),
        "steady_state" => Dict(
            "K" => s.K, "Y" => s.Y, "C" => s.C, "X" => s.X,
            "L" => s.L, "w" => s.w, "R" => s.R, "mu" => s.mu,
            "K_over_Y" => s.K_over_Y, "a_c" => s.a_c,
            "a_x" => s.a_x, "a_R" => s.a_R,
        ),
        "equilibrium_residuals" => Dict(
            String(k) => v for (k, v) in pairs(residuals)
        ),
    )

    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, data)
    end
    return path
end

function main()
    result = run_calibration()
    s = result.steady_state
    c = result.constants
    a = result.threshold
    path = normpath(joinpath(@__DIR__, "..", "results", "calibration.toml"))
    save_calibration(path, result)

    println("Nirei (2015): Issue #1 benchmark calibration")
    println("Julia version: ", VERSION)
    println("c_L = ", c.c_L, "; rho = ", c.rho)
    println("A = ", a.A, "; quadrature error estimate = ", a.integral_error)
    println("K = ", s.K, "; Y = ", s.Y, "; C = ", s.C, "; X = ", s.X)
    println("L = ", s.L, "; w = ", s.w, "; R = ", s.R, "; mu = ", s.mu)
    println("C/Y = ", s.a_c, "; X/Y = ", s.a_x)
    println("Result written: ", path)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
