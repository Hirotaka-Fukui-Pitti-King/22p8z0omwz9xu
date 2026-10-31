module Nirei2015

using Distributions: Normal, pdf, truncated
using QuadGK: quadgk
using LinearAlgebra: eigen, norm
using Random: Xoshiro, rand, rand!, randn!

export BenchmarkParameters, StructuralConstants, ThresholdConstant, SteadyState
export CalibrationResult, benchmark_parameters, lambda_distribution
export compute_structural_constants, compute_threshold_constant, compute_steady_state
export equilibrium_residuals, run_calibration

export LinearRECoefficients, LinearRESolution, ExpectedEquilibrium
export compute_linear_re_solution, compute_expected_equilibrium
export compute_saddle_path_expectations, expectation_residuals
export FirmRNGs, FirmState, firm_rngs, stream_seeds, initialize_firms
export draw_productivity!, advance_productivity!, aggregate_capital
export threshold_coefficient, realized_initial_gap_summary
export EquilibriumWorkspace, FixedPointEnumeration, ES1Selection, EquilibriumResult
export initialize_equilibrium_workspace, enumerate_fixed_points!, select_es1, solve_equilibrium!
export SimulationState, PeriodRecord, SimulationFailure, initialize_simulation, simulate_period!, simulate_path!

include("calibration.jl")
include("expectations.jl")
include("firms.jl")
include("equilibrium.jl")
include("simulation.jl")

end # module Nirei2015
