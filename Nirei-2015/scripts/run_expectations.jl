#!/usr/bin/env julia
# Run from Nirei-2015/: julia --project=. scripts/run_expectations.jl

using Nirei2015
using TOML
import Pkg

function _matrix_rows(A::Matrix{Float64})
    return [collect(A[i, :]) for i in axes(A, 1)]
end

function _dependency_versions()
    versions = Dict{String,String}()
    for info in values(Pkg.dependencies())
        if info.name in ("Distributions", "QuadGK", "TOML", "LinearAlgebra")
            versions[info.name] = isnothing(info.version) ? "stdlib" : string(info.version)
        end
    end
    return versions
end

function _expected_dict(e::ExpectedEquilibrium)
    return Dict(
        "K_prior" => e.K_prior, "mu_prior" => e.mu_prior,
        "logK_prior" => e.logK_prior, "logmu_prior" => e.logmu_prior,
        "logK_next" => e.logK_next, "logmu_next" => e.logmu_next,
        "logX_next" => e.logX_next, "logC_next" => e.logC_next,
        "logL_next" => e.logL_next, "logY_next" => e.logY_next,
        "logw_next" => e.logw_next, "logR_next" => e.logR_next,
        "K_next" => e.K_next, "mu_next" => e.mu_next,
        "X_next" => e.X_next, "C_next" => e.C_next,
        "L_next" => e.L_next, "Y_next" => e.Y_next,
        "w_next" => e.w_next, "R_next" => e.R_next,
    )
end

function main()
    calibration = run_calibration()
    solution = compute_linear_re_solution(calibration)
    steady = calibration.steady_state

    zero_state = compute_saddle_path_expectations(solution, calibration, steady.K)
    # An illustrative one-percent perturbation of *predetermined* aggregate K,
    # with the predetermined marginal-utility state on the stable eigenspace.
    perturbed_state = compute_saddle_path_expectations(
        solution, calibration, steady.K * 1.01
    )
    # A general two-state call remains available for arbitrary predetermined
    # K and mu; the simulator must not reset mu to the stable relation.
    steady_residuals = expectation_residuals(solution, calibration, zero_state)
    perturbed_residuals = expectation_residuals(solution, calibration, perturbed_state)
    residual_max = max(
        maximum(abs, Tuple(steady_residuals)),
        maximum(abs, Tuple(perturbed_residuals)),
    )
    residual_max <= 1e-10 ||
        error("Appendix A expectation equation residual exceeds 1e-10: $residual_max")

    data = Dict(
        "source" => Dict(
            "paper" => "Nirei (2015), Theoretical Economics 10, 953-985",
            "doi" => "10.3982/TE1611",
            "appendix" => "Appendix A equations (26)-(37) and Proposition 5",
            "information_timing" => "K_prior and mu_prior are observable at t=-1; outputs are E_-1 of log deviations at t=0",
            "levels" => "Exp of first-order conditional expected log deviations; not exact conditional means in levels",
            "saddle_path_example" => "The 1% K disturbance initializes mu_prior on the stable eigenspace; general state API does not reset mu",
        ),
        "environment" => Dict(
            "julia_version" => string(VERSION),
            "active_project" => something(Base.active_project(), ""),
            "dependencies" => _dependency_versions(),
        ),
        "derived_coefficients" => Dict(
            "eta_L" => solution.coefficients.eta_L,
            "A1" => solution.coefficients.A1,
            "a_R" => solution.coefficients.a_R,
        ),
        "matrices" => Dict(
            "B_rows" => _matrix_rows(solution.B),
            "D_rows" => _matrix_rows(solution.D),
            "F_rows" => _matrix_rows(solution.F),
        ),
        "eigen_solution" => Dict(
            "eigenvalues" => solution.eigenvalues,
            "stable_index" => solution.stable_index,
            "stable_eigenvalue" => solution.stable_eigenvalue,
            "unstable_eigenvalue" => solution.unstable_eigenvalue,
            "stable_eigenvector_K_normalized" => solution.stable_eigenvector,
            "h" => solution.h,
        ),
        "matrix_eigen_residuals" => Dict(
            "BF_minus_D_inf" => solution.matrix_residual,
            "stable_eigenvector_inf" => solution.eigen_residual,
            "stable_manifold_residual" => solution.saddle_residual,
        ),
        "steady_state_expectations" => _expected_dict(zero_state),
        "perturbed_expectations" => _expected_dict(perturbed_state),
        "steady_state_residuals" => Dict(string(k) => v for (k, v) in pairs(steady_residuals)),
        "perturbed_residuals" => Dict(string(k) => v for (k, v) in pairs(perturbed_residuals)),
    )

    path = normpath(joinpath(@__DIR__, "..", "results", "expectations.toml"))
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, data)
    end

    println("Nirei (2015): Issue #2 Appendix A linear RE solution")
    println("Julia version: ", VERSION)
    println("eigenvalues: ", solution.eigenvalues)
    println("stable xi: ", solution.stable_eigenvalue)
    println("saddle-path h: ", solution.h)
    println("BF-D residual: ", solution.matrix_residual)
    println("example K_{-1}: ", perturbed_state.K_prior)
    println("example E_-1[K_0] (log-linear level): ", perturbed_state.K_next)
    println("example E_-1[X_0] (log-linear level): ", perturbed_state.X_next)
    println("max expectation residual: ", residual_max)
    println("Result written: ", path)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
