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
    fixed::Bool=false,polish_households::Bool=false,
    strict_portfolio::Bool=false) where {T<:AbstractFloat}

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
                residual,candidate,ctx,w;polish_households=polish_households,
                strict_portfolio=strict_portfolio)
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
            polish_households=polish_households,
            strict_portfolio=strict_portfolio,tolerance=T(1e-6),
            maxiter=70,max_evaluations=5000,verbose=true)
    end
    return solve_stage3_flexible_full(
        ctx;x0=copy(candidate),polish_households=polish_households,
        strict_portfolio=strict_portfolio,tolerance=T(1e-6),maxiter=70,max_evaluations=5000,verbose=true)
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

"""
    baseline_guess_if_present(outputdir, label, original_ctx, fallback)

For a NEW Euler-polished inherited-portfolio experiment, a stored
baseline-EGM market root is admissible as an INITIAL AGGREGATE GUESS
only when it matches the ORIGINAL baseline boom's inherited K and B
sequence. It is never reused as a solution of the polished economy.
"""
function baseline_guess_if_present(outputdir::AbstractString,
                                   label::String,
                                   original_ctx::Stage3Context{T},
                                   fallback::Vector{T}) where {T<:AbstractFloat}
    file=joinpath(outputdir,"economic_"*label*"_MARKET_ROOT_UNCERTIFIED.csv")
    candidate=read_economic_root_guess(file,original_ctx)
    if candidate===nothing
        println("No baseline ",label," path; using computed economic guess.")
        return fallback
    end
    println("Using baseline ",label,
            " K/L/terms path ONLY as a numerical starting guess "*
            "for the re-computed Euler-polished model.")
    return candidate
end

"""
    compare_economic_baseline_if_present(outputdir, label, w)

Report changes to the first 40 economic dates against a stored baseline
market root when that root is available. Missing baseline files are
not a reason to abandon an already computed Euler-polished GE.
"""
function compare_economic_baseline_if_present(outputdir::AbstractString,
                                              label::String,w)
    file=joinpath(outputdir,"economic_"*label*"_MARKET_ROOT_UNCERTIFIED.csv")
    if !isfile(file)
        println("NO BASELINE COMPARISON for ",label,
                ": the original unpolished path was not saved.")
        return nothing
    end
    comparison=compare_boom_horizons(file,w)
    println("FIRST-40 ECONOMIC PATH RESPONSE for ",label,": ",comparison)
    return comparison
end

function main()
    valid=all(x->x in ("--boom-only","--audit-euler-refinement",
                        "--solve-polished-boom","--solve-polished-all",
                        "--stage5-solve-strict-flexible") ||
               startswith(x,"--horizon="),ARGS)
    valid || error("Usage: julia --project=. scripts/run_economic_replication.jl [--boom-only] [--horizon=300] [--audit-euler-refinement] [--solve-polished-boom] [--solve-polished-all] [--stage5-solve-strict-flexible] (choose one correction mode)")
    settings=filter(x->startswith(x,"--horizon="),ARGS)
    length(settings)<=1 || throw(ArgumentError("specify the horizon once"))
    requested_horizon=isempty(settings) ? 200 :
        parse(Int,split(only(settings),"=")[2])
    requested_horizon>=200 ||
        throw(ArgumentError("horizon cannot be below 200"))
    boom_only="--boom-only" in ARGS
    polish_audit="--audit-euler-refinement" in ARGS
    polish_solve="--solve-polished-boom" in ARGS
    polish_all="--solve-polished-all" in ARGS
    strict_solve="--stage5-solve-strict-flexible" in ARGS
    count(x->x, (polish_solve,polish_all,polish_audit,strict_solve))<=1 ||
        throw(ArgumentError("choose at most one Euler-polished run mode"))
    polish_all && (boom_only || requested_horizon!=200) &&
        throw(ArgumentError("--solve-polished-all must use the standard 200-date four-economy experiment"))
    strict_solve && (boom_only || requested_horizon!=200) &&
        throw(ArgumentError("--stage5-solve-strict-flexible requires the standard 200-date sudden-stop horizon"))
    polish_solve && !boom_only &&
        throw(ArgumentError("polished equilibrium solve requires --boom-only"))
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
    polishedbasefile=joinpath(outputdir,
        "economic_credit_expansion_euler_polished_MARKET_ROOT_UNCERTIFIED.csv")
    if requested_horizon>n0 && polish_solve
        # The 300-date corrected household model must start from its own
        # SOLVED 200-date corrected GE path. Never relabel an unpolished
        # 300-date path as a corrected-household solution.
        prior=read_economic_root_guess(polishedbasefile,originalctx)
        prior===nothing && error("Missing 200-date EULER-POLISHED market root. Run --boom-only --solve-polished-boom before extending the same model to 300 dates.")
        boomctx=FMR2020.stage3_extend_horizon(originalctx,requested_horizon)
        guess=FMR2020.stage3_prolong_aggregate_guess(
            prior,originalctx,boomctx)
        # If the user has ALREADY solved this horizon with base EGM,
        # that candidate generally lies closer to the polished root
        # than an extrapolation of the corrected 200-date equilibrium.
        # It is ONLY an initial guess: current Euler-polished market
        # residuals are independently rebuilt before accepting a root.
        unpolished_extended_csv=joinpath(outputdir,
            "economic_credit_expansion_h"*string(requested_horizon)*
            "_MARKET_ROOT_UNCERTIFIED.csv")
        extended_guess=read_economic_root_guess(
            unpolished_extended_csv,boomctx)
        if extended_guess!==nothing
            println("Using saved ",requested_horizon,
                    "-date BASE-EGM path solely as initial guess "*
                    "for the new Euler-polished GE. No base policies reused.")
            guess=extended_guess
        end
        boom_label="credit_expansion_euler_polished_h"*string(requested_horizon)
    elseif requested_horizon>n0
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
    boom=solve_or_recheck_economic(
        boomctx,guess,outputdir,boom_label;
        polish_households=polish_solve && requested_horizon>n0)
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
    if polish_solve && requested_horizon==n0
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
    if strict_solve
        # Stage 5 is a SINGLE discriminating full-GE computation:
        # original foreign-credit path and inherited pre-shock household
        # distribution are unchanged. Only the initial gross capital /
        # foreign bond decomposition is corrected to clear q_pre K and B.
        # No source MATLAB household policies or 399 equation shortcuts.
        label="credit_expansion_euler_polished"
        prex=read_economic_root_guess(polishedbasefile,boomctx)
        prex===nothing && error("Requires the archived Euler-polished boom")
        preboom=solve_or_recheck_economic(
            boomctx,prex,outputdir,label;polish_households=true)
        ctx=stage3_reversal_context(
            steady,preboom;shock_start=41,leverage_hat=1/16)
        name="flexible_leveraged_euler_polished"
        old_csv=joinpath(outputdir,
            "economic_"*name*"_MARKET_ROOT_UNCERTIFIED.csv")
        candidate=read_economic_root_guess(old_csv,ctx)
        candidate===nothing && error("Requires archived leveraged flexible-FX Euler-polished GE path")
        corrected_name="flexible_leveraged_euler_polished_STRICT_PORTFOLIO"
        println("STAGE5: resolving the 599-equation leveraged flexible "*
                "sudden stop with STRICT inherited gross capital/bonds "*
                "and identical foreign-credit shock. The old solved "*
                "path is only a NUMERICAL INITIAL GUESS.")
        newroot=solve_or_recheck_economic(
            ctx,candidate,outputdir,corrected_name;
            polish_households=true,strict_portfolio=true)
        report=economic_validation(
            newroot,ctx,corrected_name,outputdir)
        @assert report.diagnostics.max_dynamic_euler<=1e-4
        println("STAGE5 STRICT RESULT: x_pre=",ctx.terms_pre,
                ", x_impact=",newroot.workspace.terms[1],
                ", nominal_depreciation_logpoints=",
                100*log(ctx.terms_pre/newroot.workspace.terms[1]),
                ", market_residual=",
                report.diagnostics.maximum_scaled_solver_residual,
                ", Euler_max=",report.diagnostics.max_dynamic_euler,
                ", terminal_K_gap=",report.diagnostics.terminal_K_gap)
        println("A strict-portfolio root is NOT a validated alternative "*
                "equilibrium branch or a complete paper replication. "*
                "Terminal and grid convergence and an independent "*
                "second initial guess are outstanding.")
        return nothing
    end
    if polish_all
        println("Starting economically consistent four-model Euler-polished "*
                "replication; each stage is an independently solved "*
                "599-variable market system.")
        corrected_label="credit_expansion_euler_polished"
        corrected=solve_or_recheck_economic(
            boomctx,copy(boom.x),outputdir,corrected_label;
            polish_households=true)
        corrected_check=economic_validation(
            corrected,boomctx,corrected_label,outputdir)

        # Reconstruct old-context starting guesses using baseline boom
        # only for validation. The NEW shock context and portfolio
        # distribution come exclusively from corrected boom.
        old_noctx=stage3_reversal_context(
            steady,boom;shock_start=41,leverage_hat=0.0)
        new_noctx=stage3_reversal_context(
            steady,corrected;shock_start=41,leverage_hat=0.0)
        no_guess=baseline_guess_if_present(
            outputdir,"flexible_unleveraged",old_noctx,
            economic_aggregate_guess(corrected,new_noctx))
        no_corrected=solve_or_recheck_economic(
            new_noctx,no_guess,outputdir,
            "flexible_unleveraged_euler_polished";
            polish_households=true)
        no_check=economic_validation(
            no_corrected,new_noctx,
            "flexible_unleveraged_euler_polished",outputdir)
        compare_economic_baseline_if_present(
            outputdir,"flexible_unleveraged",no_corrected.workspace)

        old_levctx=stage3_reversal_context(
            steady,boom;shock_start=41,leverage_hat=1.0/16.0)
        new_levctx=stage3_reversal_context(
            steady,corrected;shock_start=41,leverage_hat=1.0/16.0)
        lev_guess=baseline_guess_if_present(
            outputdir,"flexible_leveraged",old_levctx,
            copy(no_corrected.x))
        lev_corrected=solve_or_recheck_economic(
            new_levctx,lev_guess,outputdir,
            "flexible_leveraged_euler_polished";
            polish_households=true)
        lev_check=economic_validation(
            lev_corrected,new_levctx,
            "flexible_leveraged_euler_polished",outputdir)
        compare_economic_baseline_if_present(
            outputdir,"flexible_leveraged",lev_corrected.workspace)

        sticky_terms=[new_levctx.economy.terminal.terms+
                      (new_levctx.terms_pre-new_levctx.economy.terminal.terms)*
                      0.95^(t-1) for t in 1:n]
        fixed_fallback=vcat(lev_corrected.x[1:2n-1],sticky_terms)
        fx_guess=baseline_guess_if_present(
            outputdir,"fixed_leveraged",old_levctx,fixed_fallback)
        fx_corrected=solve_or_recheck_economic(
            new_levctx,fx_guess,outputdir,
            "fixed_leveraged_euler_polished";
            fixed=true,polish_households=true)
        fx_check=economic_validation(
            fx_corrected,new_levctx,
            "fixed_leveraged_euler_polished",outputdir;fixed=true)
        compare_economic_baseline_if_present(
            outputdir,"fixed_leveraged",fx_corrected.workspace)
        println("CORRECTED FOUR-EXPERIMENT RESULTS:",
                " credit_expansion_market=",corrected_check.diagnostics.maximum_scaled_solver_residual,
                " unleveraged_market=",no_check.diagnostics.maximum_scaled_solver_residual,
                " flexible_leveraged_market=",lev_check.diagnostics.maximum_scaled_solver_residual,
                " fixed_leveraged_market=",fx_check.diagnostics.maximum_scaled_solver_residual)
        println("Results are saved as distinctly named provisional market "*
                "roots. Euler/KKT, terminal and grid checks remain binding; "*
                "Figures 1/4/5 are not certified or published.")
        return nothing
    end

    if requested_horizon>n0
        # Always compare the SAME household algorithm in both horizons.
        comparisonfile=polish_solve ? polishedbasefile : basefile
        horizon_report=compare_boom_horizons(comparisonfile,boom.workspace)
        println("SAME-SHOCK BOOM HORIZON (first 40 quarters, 200 vs ",
                n,", household_map=",
                polish_solve ? "Euler-polished" : "baseline-EGM",
                ") relative differences: ",horizon_report)
        println("Max relative early-path difference: ",
                horizon_report.maximum)
        println("TERMINAL DIAGNOSTICS for comparison: K_gap=",
                boom_check.diagnostics.terminal_K_gap,
                " terminal_budget_max=",
                boom_check.diagnostics.max_terminal_household_budget,
                " terminal_budget_mass_weighted=",
                boom_check.diagnostics.weighted_terminal_household_budget)
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
