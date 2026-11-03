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

"""
    fixed_horizon_comparison(reference_csv, longer_workspace; first_periods=40)

Compare the same early 40 quarters of an economic 200-date market root
and a longer-horizon market root, using unchanged exogenous shocks.
Reference CSVs are written by stage3_write_csv, with column names checked.
The scale for each series is max(1, |reference level|).
"""
function fixed_horizon_comparison(
    reference_csv::AbstractString, w::FMR2020.TransitionWorkspace{T};
    first_periods::Int=40) where {T<:AbstractFloat}

    isfile(reference_csv) ||
        throw(ArgumentError("missing shorter-horizon comparison: "*reference_csv))
    lines=readlines(reference_csv)
    length(lines)>first_periods ||
        throw(DimensionMismatch("shorter-horizon CSV too short"))
    names=split(strip(lines[1]),",")
    pairs=(
        ("K",w.K),("L",w.L),("Y",w.Y),
        ("w",w.wage),("r",w.rate),
        ("terms",w.terms),
        ("home_consumption",w.aggregate_home),
        ("I",w.I))
    deviations=Dict{String,T}()
    for (name,series) in pairs
        pos=findfirst(==(name),names)
        pos===nothing && throw(ArgumentError("missing CSV column "*name))
        maxchange=zero(T)
        for t in 1:first_periods
            cols=split(strip(lines[t+1]),",")
            old=parse(T,cols[pos])
            fresh=series[t]
            isfinite(old) && isfinite(fresh) ||
                throw(DomainError((old,fresh),"nonfinite horizon comparison"))
            maxchange=max(maxchange,abs(fresh-old)/max(one(T),abs(old)))
        end
        deviations[name]=maxchange
    end
    return (by_variable=deviations,maximum=maximum(values(deviations)),
            first_periods=first_periods)
end

function solve_economic_fixed_ge()
    length(ARGS) in (1,2) ||
        error("Usage: julia --project=. scripts/solve_economic_fixed_ge.jl FIXTURE_DIR [HORIZON]")
    dir=ARGS[1]
    requested_horizon=length(ARGS)==2 ? parse(Int,ARGS[2]) : 200
    requested_horizon>=200 ||
        throw(ArgumentError("author fixture starts from the 200-date experiment"))
    steady=solve_three_steady_states(default_parameters();verbose=false)
    mu=source_matrix(dir,"shock_pre_distribution")
    metadata=source_vector(dir,"shock_inherited_metadata")
    source_ctx=source_reversal_context(steady,mu,metadata,1.0/16.0)
    # A horizon extension must preserve the entire original 200-date
    # exogenous credit path; only append its stationary continuation.
    ctx=FMR2020.stage3_extend_horizon(source_ctx,requested_horizon)
    abs(ctx.theta_nr-100.0)<1e-12 ||
        error("Unexpected nominal adjustment-cost calibration")
    saved_author_x=source_vector(dir,"fixed_candidate")
    n0=length(source_ctx.economy.B)
    n=length(ctx.economy.B)
    length(saved_author_x)==3n0-1 ||
        throw(DimensionMismatch("Author reference must contain 3T0-1 saved path variables"))
    checkpoint=joinpath(dir,"julia_$(3n-1)_FIXED_TERMINAL_BC.checkpoint")
    basecheckpoint=joinpath(dir,"julia_$(3n0-1)_FIXED_TERMINAL_BC.checkpoint")
    start_x=if isfile(checkpoint)
        println("Rechecking existing Julia GE root at horizon ",n,
                " under CURRENT economic equations.")
        FMR2020.stage3_read_checkpoint(checkpoint,3n-1)
    elseif n>n0
        prior=isfile(basecheckpoint) ?
            FMR2020.stage3_read_checkpoint(basecheckpoint,3n0-1) :
            saved_author_x
        println("Extending unchanged 200-date shock to ",n,
                " dates using prior path as NUMERICAL INITIAL GUESS only.")
        FMR2020.stage3_prolong_aggregate_guess(prior,source_ctx,ctx)
    else
        println("Using original saved K/L/terms PATH as numerical guess, " *
                "never as a household policy or a certified GE.")
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
    println("FIXED ECONOMIC GE diagnostic, horizon=",n,": ",diag)
    fixed_market_acceptance(diag) ||
        error("Fixed root fails independent asset/labor/goods/price/budget/mass audit")
    provisional=joinpath(dir,"julia_$(3n-1)_FIXED_TERMINAL_BC_EULER_PENDING.csv")
    FMR2020.stage3_write_csv(result,ctx,provisional;fixed=true)
    if result.iterations>0
        FMR2020.stage3_write_checkpoint(
            checkpoint,result.x,result.residual,result.iterations)
    end
    reference_path=joinpath(dir,
        "julia_$(3n0-1)_FIXED_TERMINAL_BC_EULER_PENDING.csv")
    convergence=n>n0 && isfile(reference_path) ?
        fixed_horizon_comparison(reference_path,result.workspace) : nothing
    convergence===nothing ?
        println("HORIZON CONVERGENCE NOT YET VERIFIED. Need valid 200- and " *
                "longer-date paths with identical shock.") :
        println("SAME-SHOCK 200-to-",n,
                " DATE FIRST-40-QUARTER CONVERGENCE: ",convergence)
    horizon_pass=convergence!==nothing && convergence.maximum<=1e-3
    if !(fixed_euler_acceptance(diag) &&
         fixed_terminal_acceptance(diag,ctx.economy.terminal.K,
                                   ctx.economy.terminal.L) &&
         horizon_pass)
        println("FIXED MARKET-CLEARING ROOT SAVED: ",provisional,
                "; full dynamic GE certification pending.",
                " Euler_max=",diag.max_dynamic_euler,
                " Euler_mean=",diag.weighted_mean_dynamic_euler,
                " Euler_worst=",diag.max_dynamic_euler_state,
                " Euler_mass_above_1e-4=",diag.euler_mass_above_1e4,
                " terminal_budget_max=",diag.max_terminal_household_budget,
                " terminal_K_gap=",diag.terminal_K_gap,
                " terminal_L_gap=",diag.terminal_L_gap,
                " terminal_terms_gap=",diag.terminal_terms_gap,
                " horizon_first40_relative_gap=",
                convergence===nothing ? "not checked" : convergence.maximum)
        return nothing
    end
    certified=joinpath(dir,"julia_$(3n-1)_FIXED_TERMINAL_BC_VERIFIED.csv")
    FMR2020.stage3_write_csv(result,ctx,certified;fixed=true)
    println("FIXED GE NUMERICALLY VERIFIED, including first-40-quarter " *
            "same-shock horizon convergence: ",
            certified)
    return nothing
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    solve_economic_fixed_ge()
end
