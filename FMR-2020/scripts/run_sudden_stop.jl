
using FMR2020

function warmstart_from_stage2(filename::AbstractString,n::Int)
    isfile(filename) || return nothing
    lines=readlines(filename)
    length(lines)==n+1 || throw(ArgumentError("Stage 2 CSV horizon differs"))
    names=split(lines[1],",")
    ik=findfirst(==("K"),names); il=findfirst(==("L"),names)
    (isnothing(ik) || isnothing(il)) &&
        throw(ArgumentError("missing capital/labor in Stage 2 file"))
    K=Float64[]; L=Float64[]
    for line in lines[2:end]
        cells=split(line,",")
        push!(K,parse(Float64,cells[ik]))
        push!(L,parse(Float64,cells[il]))
    end
    return vcat(K[2:end],L)
end

function output_revaluation(ctx::Stage3Context{T},
                            flexible::FMR2020.TransitionSolveResult{T},
                            fixed::FMR2020.TransitionSolveResult{T},
                            boom::FMR2020.TransitionSolveResult{T},
                            filename::AbstractString;
                            shock_start::Int=41) where {T}
    grid=ctx.economy.grids.fine
    src=ctx.mu_pre
    K=flexible.workspace.K[1]
    rich=sum(src[i,j]*grid[i] for i in eachindex(grid),j in axes(src,2)
             if grid[i]>0)
    lev=ctx.economy.B_initial==0 ? one(T) :
         (one(T)-ctx.leverage_hat)*K*ctx.economy.q_initial/rich
    fixedcap=ctx.leverage_hat*K
    rflex=flexible.workspace.terms[1]/ctx.terms_pre
    rfixed=fixed.workspace.terms[1]/ctx.terms_pre
    # A household's inherited portfolio is the SAME in both regimes.
    # FX valuation is regime-specific, as in source backsolve_egm.m.
    function consumption_at(result,a::T,j::Int)
        if a<grid[1] || a>grid[end]
            throw(DomainError(a,"occupied revalued wealth outside grid"))
        end
        if a==grid[end]
            return result.workspace.fine_consumption[end,j,1]
        end
        k=max(1,min(length(grid)-1,searchsortedlast(grid,a)))
        frac=(a-grid[k])/(grid[k+1]-grid[k])
        left=result.workspace.fine_consumption[k,j,1]
        right=result.workspace.fine_consumption[k+1,j,1]
        # Linear interpolation is used solely for reporting household
        # consumption at their actual mapped initial wealth, NOT in
        # the nonlinear EGM or GE solver.
        if frac==zero(T)
            return left
        elseif frac==one(T)
            return right
        end
        return (one(T)-frac)*left+frac*right
    end
    open(filename,"w") do io
        println(io,"asset,productivity,mass_pre,gross_capital,foreign_bonds,asset_flexible,asset_fixed,exchange_ratio_flexible,exchange_ratio_fixed,C_before,C_flexible,C_fixed,delta_C_flexible,delta_C_fixed")
        for j in axes(src,2),i in eachindex(grid)
            mass=src[i,j]
            mass==0 && continue
            a=grid[i]
            cap=max(lev*a+fixedcap,zero(T))
            bond=a-cap
            afterflex=cap+bond/rflex
            afterfix=cap+bond/rfixed
            initial_c=boom.workspace.fine_consumption[i,j,shock_start]
            after_cflex=consumption_at(flexible,afterflex,j)
            after_cfix=consumption_at(fixed,afterfix,j)
            isfinite(initial_c) && initial_c>0 &&
              isfinite(after_cflex) && after_cflex>0 &&
              isfinite(after_cfix) && after_cfix>0 ||
                throw(DomainError((a,j,afterflex,afterfix),
                                  "invalid positive-mass revaluation consumption"))
            println(io,join((a,j,mass,cap,bond,afterflex,afterfix,
                             rflex,rfixed,initial_c,after_cflex,after_cfix,
                             after_cflex-initial_c,after_cfix-initial_c),","))
        end
    end
    return filename
end

function verify_stage3(result,ctx;fixed::Bool)
    report=stage3_diagnostics(result,ctx;fixed=fixed)
    for (metric,value) in pairs(report)
        println(fixed ? "FIXED " : "FLEX ",metric," = ",value)
    end
    good=report.max_asset_market_level<1e-4 &&
         report.max_labor_market_level<1e-4 &&
         report.max_goods_market_level<1e-4 &&
         report.max_mass_deviation<1e-8 &&
         report.max_household_budget<1e-8 &&
         report.max_nonlinear_phillips<1e-8
    good || error("Stage 3 independent GE checks failed; do not plot this candidate")
    return report
end

function main()
    println("FMR-2020 Stage 3 — foreign currency sudden stop and exchange regimes")
    println("Julia ",VERSION," / project=",Base.active_project())
    root=normpath(joinpath(@__DIR__,".."))
    out=joinpath(root,"results")
    mkpath(out)
    logpath=joinpath(out,"stage3_solver_log.txt")
    try
        steady=solve_three_steady_states(default_parameters();verbose=true)
        println("Stage 3 source audit: KK=2, foreign=true, initial k_const_hat=0")
        boomctx=stage3_expansion_context(steady)
        original=joinpath(out,"stage2_transition_expansion.csv")
        warm=warmstart_from_stage2(original,length(boomctx.economy.B))
        base_guess=isnothing(warm) ?
            initial_transition_guess(boomctx.economy) : warm
        println("STAGE3 expansion with foreign-currency denomination")
        boom=solve_stage3_flexible(boomctx;x0=base_guess,
                                   tolerance=1e-6,maxiter=60,verbose=true)
        expansion_diag=verify_stage3(boom,boomctx;fixed=false)
        reversal_ctx=stage3_reversal_context(steady,boom;
                                             shock_start=41,
                                             leverage_hat=1/16)
        inp=reversal_ctx.economy
        n=length(inp.B)
        guess=zeros(Float64,2*n-1)
        for t in 1:n-1
            guess[t]=steady.calibrated.K+
                 (inp.initial.K-steady.calibrated.K)*0.8^t
        end
        guess[n:end].=steady.calibrated.L
        println("STAGE3 surprise credit reversal — flexible exchange rate")
        flex=solve_stage3_flexible(reversal_ctx;x0=guess,
                                   tolerance=1e-6,maxiter=80,verbose=true)
        flexdiag=verify_stage3(flex,reversal_ctx;fixed=false)
        println("STAGE3 surprise credit reversal — fixed exchange rate")
        fixed=solve_stage3_fixed(reversal_ctx,flex;
                                 tolerance=1e-6,maxiter=100,verbose=true)
        fixdiag=verify_stage3(fixed,reversal_ctx;fixed=true)
        # Export ONLY the accepted economic equilibria.
        open(joinpath(out,"stage3_credit_reversal.csv"),"w") do io
            println(io,"original_quarter,B_expansion,B_after_reversal")
            for t in 1:n
                newB=t<41 ? boom.input.B[t] :
                       reversal_ctx.economy.B[t-40]
                println(io,join((t,boom.input.B[t],newB),","))
            end
        end
        stage3_write_csv(boom,boomctx,
                         joinpath(out,"stage3_corrected_expansion.csv"))
        stage3_write_csv(flex,reversal_ctx,
                         joinpath(out,"stage3_flexible_transition.csv"))
        stage3_write_csv(fixed,reversal_ctx,
                         joinpath(out,"stage3_fixed_transition.csv");fixed=true)
        output_revaluation(reversal_ctx,flex,fixed,boom,
                           joinpath(out,"stage3_household_revaluation.csv"))
        diagfile=joinpath(out,"stage3_equilibrium_diagnostics.csv")
        open(diagfile,"w") do io
            println(io,"scenario,metric,value")
            for (name,d) in (("expansion",expansion_diag),
                             ("flexible",flexdiag),("fixed",fixdiag))
                for (k,v) in pairs(d)
                    println(io,join((name,k,v),","))
                end
            end
        end
        f1=stage3_figure_1(boom,reversal_ctx.economy.B,
                           joinpath(out,"figure_1.svg"))
        f45=stage3_figures_4_5(boom,flex,fixed,out)
        open(logpath,"w") do io
            println(io,"Stage3 verification on Julia ",VERSION)
            for (name,sol,di) in (("expansion",boom,expansion_diag),
                                  ("flexible",flex,flexdiag),
                                  ("fixed",fixed,fixdiag))
                println(io,name,": converged=",sol.converged,
                           " max=",sol.maximum_scaled_residual,
                           " iterations=",sol.iterations,
                           " evaluations=",sol.evaluations)
                for (k,v) in pairs(di)
                    println(io,string(k),"=",v)
                end
            end
            println(io,"Source-based plots: ",f1,", ",join(f45,", "))
            println(io,"Actual author numerical-series match: NOT VERIFIED")
        end
        println("STAGE3 ACCEPTED: independent checks passed; files in ",out)
    catch err
        open(logpath,"w") do io
            println(io,"STAGE3 FAILED: ",sprint(showerror,err))
            println(io,"No solution or figure may be claimed from this failure.")
        end
        rethrow()
    end
end

main()
