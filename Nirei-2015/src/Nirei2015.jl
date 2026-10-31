module Nirei2015

using Distributions: Normal, pdf, truncated
using QuadGK: quadgk
using LinearAlgebra: eigen, norm

export BenchmarkParameters, StructuralConstants, ThresholdConstant, SteadyState
export CalibrationResult, benchmark_parameters, lambda_distribution
export compute_structural_constants, compute_threshold_constant, compute_steady_state
export equilibrium_residuals, run_calibration

export LinearRECoefficients, LinearRESolution, ExpectedEquilibrium
export compute_linear_re_solution, compute_expected_equilibrium
export compute_saddle_path_expectations, expectation_residuals

include("calibration.jl")
include("expectations.jl")

end # module Nirei2015
