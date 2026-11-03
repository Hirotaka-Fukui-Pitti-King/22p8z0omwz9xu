# Original economic model (NOT the source MATLAB 399/599 residual shortcuts):
# a solved boom generates the unexpected-shock inherited distribution;
# an independent full goods-clearing GE is then solved in both regimes.
# This is an expensive user-PC research calculation, NEVER run by CI.
#
# From FMR-2020 (Windows PowerShell):
#   julia --project=. scripts\run_economic_replication.jl
#
# The first 40 periods of market-clearing roots are NOT plotted as verified
# Figures 1/4/5 unless the household Euler/KKT and terminal checks pass.

using FMR2020

function economic_aggregate_guess(
    boom::FMR2020.TransitionSolveResult{T},
    ctx::Stage3Context{T};shock_start::Int=41) where {T<:AbstractFloat}
    n=length(ctx.economy.B)
    nb=length(boom.input.B)
    2<=shock_start<=nb ||
        throw(ArgumentError("surprise shock must have an existing pre-path"))
    K=Vector{T}(undef,n)
    L=similar(K)
    terms=similar(K)
    @inbounds for t in 1:n
        oldt=shock_start+t-1
        if oldt<=nb
            K[t]=boom.workspace.K[oldt]
            L[t]=boom.workspace.L[oldt]
            terms[t]=boom.workspace.terms[oldt]
        else
            K[t]=ctx.economy.terminal.K
            L[t]=ctx.economy.terminal.L
            terms[t]=ctx.economy.terminal.terms
        end
    end
    abs(K[1]-ctx.economy.initial.K)<=T(1e-8) ||
        throw(ArgumentError("shock date capital is not inherited from boom"))
    return vcat(K[2:end],L,terms)
end

"""
Recover K[2:T], L[1:T], x[1:T] from a saved JULIA market-root CSV.
The credit sequence, dates and inherited K are verified before reuse.
Only the numerical aggregate starting path is retained.
"""
function read_economic_root_guess(path::AbstractString,
                                  ctx::Stage3Context{T}) where {T<:AbstractFloat}
    isfile(path) || return nothing
    lines=readlines(path)
    n=length(ctx.economy.B)
    length(lines)==n+1 ||
        throw(DimensionMismatch("saved GE path horizon differs"))
    columns=strip.(split(lines[1],","))
    indices=map(("period","B","K","L","terms")) do key
        idx=findfirst(==(key),columns)
        isnothing(idx) && throw(ArgumentError("missing saved GE column "*key))
        idx
    end
    K=Vector{T}(undef,n)
    L=similar(K)
    x=similar(K)
    for t in 1:n
        values=split(lines[t+1],",")
        parse(Int,values[indices[1]])==t ||
            throw(ArgumentError("saved GE date mismatch"))
        b=parse(T,values[indices[2]])
        isapprox(b,ctx.economy.B[t];rtol=T(1e-10),atol=T(1e-10)) ||
            throw(ArgumentError("saved GE uses a different foreign-credit shock"))
        K[t]=parse(T,values[indices[3]])
        L[t]=parse(T,values[indices[4]])
        x[t]=parse(T,values[indices[5]])
    end
    isapprox(K[1],ctx.economy.initial.K;atol=T(1e-10),rtol=T(1e-10)) ||
        throw(ArgumentError("saved GE uses different inherited capital"))
    z=vcat(K[2:end],L,x)
    all(v->isfinite(v)&&v>zero(T),z) ||
        throw(DomainError(z,"nonpositive saved aggregate path"))
    return z
end

function economic_validation(result::FMR2020.TransitionSolveResult{T},
                             ctx::Stage3Context{T},
                             name::String,outputdir::String;
                             fixed::Bool=false) where {T}
    report=stage3_diagnostics(result,ctx;fixed=fixed,
                              explicit_flexible=!fixed)
    for (key,value) in pairs(report)
        println(name," ",key,"=",value)
    end
    market=(report.maximum_scaled_solver_residual<=T(1e-6) &&
            report.max_asset_market_level<=T(1e-4) &&
            report.max_labor_market_level<=T(1e-4) &&
            report.max_goods_market_level<=T(1e-4) &&
            report.max_household_budget<=T(1e-8) &&
            report.max_mass_deviation<=T(1e-8) &&
            report.max_nonlinear_phillips<=T(1e-8))
    market || error("Economic "*name*" market/budget/mass verification FAILED")
    # Preserve useful roots, but never mislabel them as infinite-horizon GE.
    path=joinpath(outputdir,"economic_"*name*"_MARKET_ROOT_UNCERTIFIED.csv")
    stage3_write_csv(result,ctx,path;fixed=fixed)
    println("PRESERVED UNVERIFIED ECONOMIC MARKET ROOT: ",path)
    terminal=(report.max_terminal_household_budget<=T(1e-6) &&
              report.terminal_K_gap/
                  max(one(T),abs(ctx.economy.terminal.K))<=T(1e-4) &&
              report.terminal_L_gap/
                  max(one(T),abs(ctx.economy.terminal.L))<=T(1e-4) &&
              report.terminal_terms_gap<=T(1e-6))
    euler=(report.max_dynamic_euler<=T(1e-4) &&
           report.weighted_mean_dynamic_euler<=T(1e-5))
    if !(euler&&terminal)
        println("UNVERIFIED ",name,": market root only. "*
                "Euler/KKT and terminal checks remain unpassed. "*
                "Continue shock experiments provisionally; suppress final figures.")
    else
        println("ONE-HORIZON PASS ",name,
                "; grid and horizon convergence still unverified.")
    end
    return (diagnostics=report,one_horizon_pass=euler&&terminal)
end

function main()
    isempty(ARGS) ||
        error("Usage: julia --project=. scripts/run_economic_replication.jl")
    base=normpath(joinpath(@__DIR__,".."))
    outputdir=joinpath(base,"results","economic_reproduction")
    mkpath(outputdir)
    steady=solve_three_steady_states(default_parameters();verbose=false)
    boomctx=stage3_expansion_context(steady)
    n=length(boomctx.economy.B)
    guess=vcat(initial_transition_guess(boomctx.economy),
               fill(steady.initial.terms,n))
    println("Solving corrected full-goods 599-unknown initial credit expansion")
    boom=solve_stage3_flexible_full(
        boomctx;x0=guess,tolerance=1e-6,maxiter=70,
        max_evaluations=5000,verbose=true)
    boom_check=economic_validation(boom,boomctx,"credit_expansion",outputdir)

    # The PRE-SHOCK wealth distribution, K, q and terms all come from the
    # computed Julia boom; they are never copied from archived MATLAB MAT.
    noctx=stage3_reversal_context(steady,boom;
                                  shock_start=41,leverage_hat=0.0)
    noguess=economic_aggregate_guess(boom,noctx)
    println("Solving unexpected credit reversal without household leverage")
    nolev=solve_stage3_flexible_full(
        noctx;x0=noguess,tolerance=1e-6,maxiter=70,
        max_evaluations=5000,verbose=true)
    nolev_check=economic_validation(nolev,noctx,"flexible_unleveraged",outputdir)

    ctx=stage3_reversal_context(steady,boom;
                                shock_start=41,leverage_hat=1.0/16.0)
    println("Solving foreign-currency leveraged reversal, flexible prices")
    flexible=solve_stage3_flexible_full(
        ctx;x0=copy(nolev.x),tolerance=1e-6,maxiter=70,
        max_evaluations=5000,verbose=true)
    flexible_check=economic_validation(flexible,ctx,"flexible_leveraged",outputdir)

    # Fixed exchange rates: same current shocks and inherited balance
    # sheets, but theta_nr=100 Rotemberg costs and monetary rule applied.
    # The starting guess is constructed entirely from the computed
    # flexible-price path and the inherited pre-shock terms.
    fxguess=vcat(flexible.x[1:2n-1],
                 fill(ctx.terms_pre,n))
    println("Solving fixed nominal exchange rate with nonlinear Rotemberg")
    fixed=FMR2020.stage3_nonlinear_solve(
        ctx,fxguess;fixed=true,tolerance=1e-6,maxiter=70,
        max_evaluations=5000,verbose=true)
    fixed_check=economic_validation(fixed,ctx,"fixed_leveraged",outputdir;
                        fixed=true)

    # This output is explicitly marked as PROVISIONAL until independent
    # 200-vs-300 horizon and asset-grid stability checks are completed.
    all_pass=all((boom_check.one_horizon_pass,
                  nolev_check.one_horizon_pass,
                  flexible_check.one_horizon_pass,
                  fixed_check.one_horizon_pass))
    if !all_pass
        println("FOUR MARKET ROOTS SAVED, BUT NOT VERIFIED "*
                "AS INFINITE-HORIZON EQUILIBRIA. Figures suppressed.")
        return nothing
    end
    provisionaldir=joinpath(outputdir,"PROVISIONAL_NOT_HORIZON_GRID_VALIDATED")
    mkpath(provisionaldir)
    figures=(
        stage3_figure_1(boom,ctx.economy.B,
                        joinpath(provisionaldir,"figure_1.svg")),
        stage3_figures_4_5(boom,flexible,fixed,provisionaldir)...)
    println("PROVISIONAL FIGURES GENERATED (200-date candidates only): ",
            figures)
    println("Figures are NOT a completed replication until the same-shock "*
            "horizon and asset-grid convergence checks pass.")
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end
