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
            # The initial guess must not jump at the last available
            # credit-expansion observation. Geometrically splice toward
            # the stationary tail; this changes guesses, never equations.
            tail=oldt-nb
            weight=T(0.9)^tail
            K[t]=ctx.economy.terminal.K+
                (boom.workspace.K[end]-ctx.economy.terminal.K)*weight
            L[t]=ctx.economy.terminal.L+
                (boom.workspace.L[end]-ctx.economy.terminal.L)*weight
            terms[t]=ctx.economy.terminal.terms+
                (boom.workspace.terms[end]-ctx.economy.terminal.terms)*weight
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

"""
Recompute the full 599 residual at a saved Julia path before reuse.
When the current economic equations are different, solve afresh from
that candidate rather than assuming a stored convergence flag.
"""
function solve_or_recheck_economic(
    ctx::Stage3Context{T}, initial_x::Vector{T},
    directory::AbstractString, label::String;
    fixed::Bool=false,polish_households::Bool=false) where {T<:AbstractFloat}

    n=length(ctx.economy.B)
    length(initial_x)==3n-1 ||
        throw(DimensionMismatch("full GE requires 3T-1 aggregate unknowns"))
    saved=joinpath(directory,"economic_"*label*"_MARKET_ROOT_UNCERTIFIED.csv")
    earlier=read_economic_root_guess(saved,ctx)
    candidate=earlier===nothing ? initial_x : earlier
    if earlier!==nothing
        w=FMR2020.TransitionWorkspace(ctx.economy)
        residual=zeros(T,3n-1)
        if fixed
            FMR2020.stage3_fixed_residual!(
                residual,candidate,ctx,w;polish_households=polish_households)
        else
            FMR2020.stage3_flexible_full_residual!(
                residual,candidate,ctx,w;polish_households=polish_households)
        end
        maxerror=maximum(abs,residual)
        println("RECHECK SAVED ",label,
                " CURRENT 599-D MARKET MAX=",maxerror)
        if maxerror<=T(1e-6)
            return FMR2020.TransitionSolveResult(
                ctx.economy,w,copy(candidate),residual,0,
                w.function_evals,maxerror,true)
        end
        println("Old numerical candidate does not solve CURRENT equations.")
    end
    if fixed
        return FMR2020.stage3_nonlinear_solve(
            ctx,copy(candidate);fixed=true,
            polish_households=polish_households,tolerance=T(1e-6),
            maxiter=70,max_evaluations=5000,verbose=true)
    end
    return solve_stage3_flexible_full(
        ctx;x0=copy(candidate),polish_households=polish_households,
        tolerance=T(1e-6),maxiter=70,max_evaluations=5000,verbose=true)
end

function economic_validation(result::FMR2020.TransitionSolveResult{T},
                             ctx::Stage3Context{T},
                             name::String,outputdir::String;
                             fixed::Bool=false) where {T}
    # Save a numerically converged market-root CANDIDATE before applying
    # stricter household/terminal diagnostics. An independent audit can
    # throw (e.g., off-grid optimal saving); this must not destroy an
    # expensive computed candidate or present it as a certified GE.
    result.converged || throw(ArgumentError("unconverged GE cannot be exported"))
    path=joinpath(outputdir,"economic_"*name*"_MARKET_ROOT_UNCERTIFIED.csv")
    stage3_write_csv(result,ctx,path;fixed=fixed)
    println("CANDIDATE SAVED, economic verification still pending: ",path)
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
    println("MARKET AUDIT PASSED for provisional root: ",path)
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

"""
Compare the SAME 40 shock/credit-expansion dates in a saved 200-period
Julia market root and the genuinely recomputed extended-horizon root.
The stationary extension keeps the initial shock EXACTLY unchanged.
Returned gaps are max |extended - 200-period|/max(1,|200-period|).
"""
function compare_boom_horizons(reference_csv::AbstractString,w;
                               first_periods::Int=40)
    first_periods>=1 ||
        throw(ArgumentError("early-path comparison must be nonempty"))
    isfile(reference_csv) ||
        throw(ArgumentError("missing 200-date benchmark market root"))
    lines=readlines(reference_csv)
    length(lines)>first_periods ||
        throw(DimensionMismatch("200-date benchmark is too short"))
    headers=strip.(split(lines[1],","))
    series=(("K",w.K),("L",w.L),("Y",w.Y),("I",w.I),
            ("w",w.wage),("r",w.rate),
            ("terms",w.terms),("home_consumption",w.aggregate_home))
    maxbyname=Dict{String,Float64}()
    for (name,path) in series
        length(path)>=first_periods ||
            throw(DimensionMismatch("extended path is too short"))
        index=findfirst(==(name),headers)
        isnothing(index) && throw(ArgumentError("missing benchmark column "*name))
        gap=0.0
        for t in 1:first_periods
            cells=strip.(split(lines[t+1],","))
            baseline=parse(Float64,cells[index])
            fresh=Float64(path[t])
            isfinite(baseline) && isfinite(fresh) ||
                throw(DomainError((baseline,fresh),"invalid first-40 path"))
            gap=max(gap,abs(fresh-baseline)/max(1.0,abs(baseline)))
        end
        maxbyname[name]=gap
    end
    return (maximum=maximum(values(maxbyname)),
            by_variable=maxbyname,first_periods=first_periods)
end

function main()
    valid=all(x->x in ("--boom-only","--audit-euler-refinement",
                        "--solve-polished-boom") ||
               startswith(x,"--horizon="),ARGS)
    valid || error("Usage: julia --project=. scripts/run_economic_replication.jl [--boom-only] [--horizon=300] [--audit-euler-refinement] [--solve-polished-boom]")
    settings=filter(x->startswith(x,"--horizon="),ARGS)
    length(settings)<=1 || throw(ArgumentError("specify the horizon once"))
    requested_horizon=isempty(settings) ? 200 :
        parse(Int,split(only(settings),"=")[2])
    requested_horizon>=200 ||
        throw(ArgumentError("horizon cannot be below 200"))
    boom_only="--boom-only" in ARGS
    polish_audit="--audit-euler-refinement" in ARGS
    polish_solve="--solve-polished-boom" in ARGS
    polish_audit && polish_solve &&
        throw(ArgumentError("choose either audit or polished equilibrium solve"))
    polish_solve && (!boom_only || requested_horizon!=200) &&
        throw(ArgumentError("polished equilibrium solve requires --boom-only and 200 dates"))
    polish_audit && (!boom_only || requested_horizon!=200) &&
        throw(ArgumentError("Euler refinement audit requires --boom-only and the 200-date economic root"))
    requested_horizon>200 && !boom_only &&
        throw(ArgumentError("extended-horizon calculation must be --boom-only to prevent silently overwriting the 200-date surprise experiments"))
    base=normpath(joinpath(@__DIR__,".."))
    outputdir=joinpath(base,"results","economic_reproduction")
    mkpath(outputdir)
    steady=solve_three_steady_states(default_parameters();verbose=false)
    originalctx=stage3_expansion_context(steady)
    n0=length(originalctx.economy.B)
    baseguess=vcat(initial_transition_guess(originalctx.economy),
                   fill(steady.initial.terms,n0))
    basefile=joinpath(outputdir,
        "economic_credit_expansion_MARKET_ROOT_UNCERTIFIED.csv")
    if requested_horizon>n0
        prior=read_economic_root_guess(basefile,originalctx)
        prior===nothing && error("Missing 200-date Julia boom market root. Run the default economic replication once before requesting --horizon=300.")
        boomctx=FMR2020.stage3_extend_horizon(originalctx,requested_horizon)
        guess=FMR2020.stage3_prolong_aggregate_guess(
            prior,originalctx,boomctx)
        boom_label="credit_expansion_h"*string(requested_horizon)
    else
        boomctx=originalctx
        guess=baseguess
        boom_label="credit_expansion"
    end
    n=length(boomctx.economy.B)
    println("Solving full goods-clearing initial credit expansion T=",n)
    boom=solve_or_recheck_economic(boomctx,guess,outputdir,boom_label)
    boom_check=economic_validation(boom,boomctx,boom_label,outputdir)
    if polish_audit
        # A separate workspace protects the saved market root. Perform
        # exactly one fixed-aggregate-policy correction, forward transport,
        # and market reevaluation. NO false GE-convergence claim.
        scratch=FMR2020.TransitionWorkspace(boomctx.economy)
        old_res=zeros(Float64,3n-1)
        FMR2020.stage3_flexible_full_residual!(
            old_res,boom.x,boomctx,scratch)
        improvement=FMR2020.stage3_euler_polish_households!(
            scratch,boomctx)
        FMR2020.stage3_distribution_forward!(scratch,boomctx)
        refined_res=zeros(Float64,3n-1)
        FMR2020.stage3_flexible_full_market_residual!(
            refined_res,scratch,boomctx)
        println("EULER POLICY REFINEMENT (fixed aggregate prices): ",
                improvement)
        println("BEFORE POLISH MAX SCALED 599 RESIDUAL: ",
                maximum(abs,old_res))
        println("AFTER POLISH MAX SCALED 599 RESIDUAL: ",
                maximum(abs,refined_res))
        println("MAX MARKET RESIDUAL SHIFT: ",
                maximum(abs,refined_res.-old_res))
        println("This is a household-accuracy experiment, NOT an accepted "*
                "new aggregate equilibrium. The coarse/fine policies and "*
                "aggregate markets must be solved together if the shift "*
                "is non-negligible.")
        return nothing
    end
    if polish_solve
        # Compute NEW GE root with the corrected household map INSIDE
        # every market residual and every finite-difference Jacobian.
        # The old unpolished root is a starting guess, not an accepted
        # market-clearing result under the corrected household operator.
        corrected_label="credit_expansion_euler_polished"
        old_x=copy(boom.x)
        corrected=solve_or_recheck_economic(
            boomctx,old_x,outputdir,corrected_label;
            polish_households=true)
        revised=economic_validation(
            corrected,boomctx,corrected_label,outputdir)
        println("POLISHED GE MARKET MAX SCALED RESIDUAL=",
                revised.diagnostics.maximum_scaled_solver_residual)
        println("POLISHED GE EULER MAX=",
                revised.diagnostics.max_dynamic_euler,
                " first40=",revised.diagnostics.max_dynamic_euler_first40,
                " last20=",revised.diagnostics.max_dynamic_euler_last20)
        println("POLISHED GE TERMINAL BUDGET MAX=",
                revised.diagnostics.max_terminal_household_budget,
                " K_GAP=",revised.diagnostics.terminal_K_gap)
        # Compare actual economic early paths under the unpolished and
        # polished policy operators at the SAME 200-date shock.
        change=compare_boom_horizons(basefile,corrected.workspace)
        println("FIRST-40 ECONOMIC PATH RESPONSE TO HOUSEHOLD EULER "*
                "CORRECTION (same T, same shock): ",change)
        println("POLISHED GE PATH SAVED SEPARATELY, "*
                "grid/horizon and terminal checks remain independent. "*
                "Original unpolished market paths are unchanged.")
        return nothing
    end
    if requested_horizon>n0
        horizon_report=compare_boom_horizons(basefile,boom.workspace)
        println("SAME-SHOCK BOOM HORIZON (first 40 quarters, 200 vs ",
                n,") relative differences: ",horizon_report)
        println("Max relative early-path difference: ",
                horizon_report.maximum,
                " (horizon approximation requires an additional refinement if not small)")
    end
    if boom_only
        println("BOOM-ONLY economic market-root experiment finished. "*
                "No reversal or final figures attempted.")
        return nothing
    end

    # The PRE-SHOCK wealth distribution, K, q and terms all come from the
    # computed Julia boom; they are never copied from archived MATLAB MAT.
    noctx=stage3_reversal_context(steady,boom;
                                  shock_start=41,leverage_hat=0.0)
    noguess=economic_aggregate_guess(boom,noctx)
    println("Solving unexpected credit reversal without household leverage")
    nolev=solve_or_recheck_economic(
        noctx,noguess,outputdir,"flexible_unleveraged")
    nolev_check=economic_validation(nolev,noctx,"flexible_unleveraged",outputdir)

    ctx=stage3_reversal_context(steady,boom;
                                shock_start=41,leverage_hat=1.0/16.0)
    println("Solving foreign-currency leveraged reversal, flexible prices")
    flexible=solve_or_recheck_economic(
        ctx,copy(nolev.x),outputdir,"flexible_leveraged")
    flexible_check=economic_validation(flexible,ctx,"flexible_leveraged",outputdir)

    # Fixed exchange rates: same current shocks and inherited balance
    # sheets, but theta_nr=100 Rotemberg costs and monetary rule applied.
    # The starting guess is constructed entirely from the computed
    # flexible-price path and the inherited pre-shock terms.
    # Rotemberg pricing can reject the flexible-ToT surprise jump:
    # theta_nr=100 would imply extreme or negative real marginal cost.
    # Start with a slow nominal-price adjustment from inherited ToT
    # while retaining already solved flexible K and L as numerical seeds.
    sticky_terms=[ctx.economy.terminal.terms+
                  (ctx.terms_pre-ctx.economy.terminal.terms)*
                  0.95^(t-1) for t in 1:n]
    fxguess=vcat(flexible.x[1:2n-1],sticky_terms)
    println("Solving fixed nominal exchange rate with nonlinear Rotemberg")
    fixed=solve_or_recheck_economic(
        ctx,fxguess,outputdir,"fixed_leveraged";fixed=true)
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
