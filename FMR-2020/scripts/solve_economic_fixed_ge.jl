# Solve the ECONOMIC (not source-implementation) fixed-exchange-rate GE.
# Original author MATLAB solution is an INITIAL GUESS ONLY; its numerical
# savings policies, interpolation, consumption floor and nonzero residuals
# are NEVER copied into Julia's equilibrium equations.
#
# From FMR-2020 on the local PC:
# julia --project=. scripts/solve_economic_fixed_ge.jl results/author_reference
#
# No claims of complete GE unless both market and independent Euler/KKT
# conditions pass. A market-clearing path with a fine-grid Euler error is
# preserved separately, with an explicit provisional label.

using FMR2020
using LinearAlgebra
include(joinpath(@__DIR__,"evaluate_author_explicit_goods.jl"))

function fixed_market_acceptance(diag)
    return diag.maximum_scaled_solver_residual<=1e-6 &&
           diag.max_asset_market_level<=1e-4 &&
           diag.max_labor_market_level<=1e-4 &&
           diag.max_goods_market_level<=1e-4 &&
           diag.max_household_budget<=1e-8 &&
           diag.max_mass_deviation<=1e-8 &&
           diag.max_nonlinear_phillips<=1e-8
end

function fixed_euler_acceptance(diag)
    return diag.max_dynamic_euler<=1e-4 &&
           diag.weighted_mean_dynamic_euler<=1e-5
end

function fixed_terminal_acceptance(diag,terminal_K::Float64,terminal_L::Float64)
    return diag.max_terminal_household_budget<=1e-6 &&
           diag.terminal_K_gap/max(1.0,abs(terminal_K))<=1e-4 &&
           diag.terminal_L_gap/max(1.0,abs(terminal_L))<=1e-4 &&
           diag.terminal_terms_gap<=1e-6
end

function solve_economic_fixed_ge()
    length(ARGS)==1 ||
        error("Usage: julia --project=. scripts/solve_economic_fixed_ge.jl FIXTURE_DIR")
    dir=ARGS[1]
    steady=solve_three_steady_states(default_parameters();verbose=false)
    mu=source_matrix(dir,"shock_pre_distribution")
    metadata=source_vector(dir,"shock_inherited_metadata")
    ctx=source_reversal_context(steady,mu,metadata,1.0/16.0)
    # theta_nr is original economic calibration; all policy equations,
    # borrowing feasibility, and distribution transport are Julia's own.
    abs(ctx.theta_nr-100.0)<1e-12 ||
        error("Unexpected nominal adjustment-cost calibration")
    saved_author_x=source_vector(dir,"fixed_candidate")
    n=length(ctx.economy.B)
    length(saved_author_x)==3n-1 ||
        throw(DimensionMismatch("Fixed exchange GE needs 3T-1 path variables"))
    checkpoint=joinpath(dir,"julia_599_FIXED_TERMINAL_BC.checkpoint")
    start_x=if isfile(checkpoint)
        println("Using previously saved Julia fixed-market solution as guess; " *
                "all economic residuals will be evaluated again.")
        FMR2020.stage3_read_checkpoint(checkpoint,3n-1)
    else
        println("Using MATLAB K/L/terms PATH ONLY as starting guess, " *
                "NOT MATLAB policies or MATLAB residuals.")
        saved_author_x
    end
    w=FMR2020.TransitionWorkspace(ctx.economy)
    initial_res=zeros(Float64,3n-1)
    FMR2020.stage3_fixed_residual!(initial_res,start_x,ctx,w)
    println("FIXED ECONOMIC GE initial scaled max residual=",
            maximum(abs,initial_res),
            " initial L2=",norm(initial_res))
    result=FMR2020.stage3_nonlinear_solve(
        ctx,copy(start_x);fixed=true,tolerance=1e-6,maxiter=70,
        max_evaluations=5000,verbose=true)
    result.converged ||
        error("Fixed exchange market residual did not converge")
    diag=FMR2020.stage3_diagnostics(result,ctx;fixed=true)
    println("FIXED ECONOMIC GE diagnostic: ",diag)
    fixed_market_acceptance(diag) ||
        error("Fixed root fails independent asset/labor/goods/price/budget/mass audit")
    provisional=joinpath(dir,"julia_599_FIXED_TERMINAL_BC_EULER_PENDING.csv")
    FMR2020.stage3_write_csv(result,ctx,provisional;fixed=true)
    if result.iterations>0
        FMR2020.stage3_write_checkpoint(
            checkpoint,result.x,result.residual,result.iterations)
    end
    if !(fixed_euler_acceptance(diag) &&
         fixed_terminal_acceptance(diag,ctx.economy.terminal.K,
                                   ctx.economy.terminal.L))
        println("FIXED MARKET-CLEARING ROOT SAVED: ",provisional,
                "; full dynamic GE certification pending.",
                " Euler_max=",diag.max_dynamic_euler,
                " Euler_mean=",diag.weighted_mean_dynamic_euler,
                " Euler_worst=",diag.max_dynamic_euler_state,
                " Euler_mass_above_1e-4=",diag.euler_mass_above_1e4,
                " terminal_budget_max=",diag.max_terminal_household_budget,
                " terminal_K_gap=",diag.terminal_K_gap,
                " terminal_L_gap=",diag.terminal_L_gap,
                " terminal_terms_gap=",diag.terminal_terms_gap)
        return nothing
    end
    certified=joinpath(dir,"julia_599_FIXED_TERMINAL_BC_VERIFIED.csv")
    FMR2020.stage3_write_csv(result,ctx,certified;fixed=true)
    println("FIXED FULL GE VERIFIED (market, budget, mass, NKPC, Euler): ",
            certified)
    return nothing
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    solve_economic_fixed_ge()
end
