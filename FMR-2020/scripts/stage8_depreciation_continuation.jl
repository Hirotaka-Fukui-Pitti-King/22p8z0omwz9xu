# Stage 8 — one source-informed candidate in the EXISTING strict-portfolio GE.
#
# Run from FMR-2020:
#   julia --project=. scripts/stage8_depreciation_continuation.jl
#
# The default is ONE 599-vector residual evaluation, with at most three
# perturbed evaluations if the source candidate is feasible. No GE solve.
# A bounded full nonlinear attempt is opt-in, on the owner's Windows PC:
#   julia --project=. scripts/stage8_depreciation_continuation.jl --solve
#
# The original source 399-vector is ONLY a numerical starting guess:
# the inherited distribution, capital, B path, terminal policy, and
# strict gross capital/foreign bond accounting are from the SAME
# archived, Euler-polished Julia boom that generated the Stage 5 root.
# No original MATLAB source workspace or previously solved root is edited.

using FMR2020
using LinearAlgebra

include(joinpath(@__DIR__,"stage5_initial_fx_audit.jl"))
include(joinpath(@__DIR__,"evaluate_author_explicit_goods.jl"))

"""
    stage8_author_initial_guess(reference_dir, ctx)

Numerical initialization only. Reuse the saved original MATLAB leveraged
399-vector plus 200 terms, but inherit K1 and the entire economy from ctx.
"""
function stage8_author_initial_guess(reference_dir::AbstractString,ctx::Stage3Context{T}) where {T<:AbstractFloat}
    n=length(ctx.economy.B)
    original=source_vector(reference_dir,"leveraged_candidate")
    prices=source_vector(reference_dir,"leveraged_source_terms")
    meta=source_vector(reference_dir,"shock_inherited_metadata")
    n==200 && length(original)==2n-1 && length(prices)==n ||
        throw(DimensionMismatch("Stage 8 requires original 399/200 source benchmark"))
    length(meta)==7 || throw(DimensionMismatch("source inherited metadata"))
    x=T.(vcat(original,prices))
    all(v->isfinite(v)&&v>zero(T),x) ||
        throw(DomainError(x,"source-initialized aggregate guess must be positive"))
    return x,meta
end

"""
    stage8_portfolio_support(prices_work, ctx, x1)

Check financial viability of actual positive-mass PRE-revaluation states,
before lottery interpolation. Uses the same strict gross portfolio scale
as the full GE, but never changes the source candidate or the model.
"""
function stage8_portfolio_support(w::FMR2020.TransitionWorkspace{T},
                                  ctx::Stage3Context{T},x1::T) where {T<:AbstractFloat}
    inp=ctx.economy
    mu=ctx.mu_pre
    grid=inp.grids.fine
    priorgross=zero(T)
    poswealth=zero(T)
    oldwealth=zero(T)
    for j in axes(mu,2),i in eachindex(grid)
        mass=mu[i,j]
        mass==zero(T) && continue
        oldwealth+=mass*grid[i]
        if grid[i]>zero(T)
            poswealth+=mass*grid[i]
        end
    end
    lev=inp.B_initial==zero(T) ? one(T) :
        (one(T)-ctx.leverage_hat)*inp.initial.K*inp.q_initial/poswealth
    capital_offset=ctx.leverage_hat*inp.initial.K
    for j in axes(mu,2),i in eachindex(grid)
        mass=mu[i,j]
        mass==zero(T) && continue
        priorgross+=mass*max(lev*grid[i]+capital_offset,zero(T))
    end
    target=inp.q_initial*inp.initial.K
    scale=target/priorgross
    ratio=ctx.terms_pre/x1
    badmass=zero(T)
    badstates=0
    minresource=T(Inf)
    minstate=(0,0)
    minwealth=T(Inf)
    outgridmass=zero(T)
    for j in axes(mu,2),i in eachindex(grid)
        mass=mu[i,j]
        mass==zero(T) && continue
        grosscapital=scale*max(lev*grid[i]+capital_offset,zero(T))
        foreignbonds=grid[i]-grosscapital
        wealth=grosscapital+ratio*foreignbonds
        earnings=w.wage[1]*w.L[1]*inp.chain.levels[j]+w.profits[1]
        spendable=(one(T)+w.rate[1])*wealth+earnings-inp.grids.borrowing_limit
        if spendable<minresource
            minresource=spendable
            minstate=(i,j)
            minwealth=wealth
        end
        if spendable<=zero(T)
            badmass+=mass
            badstates+=1
        end
        if wealth<grid[1]-T(1e-9) || wealth>grid[end]+T(1e-9)
            outgridmass+=mass
        end
    end
    return (minimum_resources=minresource,
            worst_asset_index=minstate[1],
            worst_income_index=minstate[2],
            worst_revalued_wealth=minwealth,
            violating_mass=badmass,
            violating_states=badstates,
            offgrid_mass=outgridmass,
            gross_capital_gap=priorgross-target,
            strict_gross_scale=scale,
            inherited_wealth_gap=oldwealth-target-inp.B_initial)
end

function stage8_store_metrics(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"metric,value")
        for (key,val) in rows
            value=replace(string(val),","=>";", "\n"=>" ")
            println(io,key,",",value)
        end
    end
    return path
end

"""
    stage8_independent_micro_stats(work, ctx)

Verify the independently reconstructed positive-mass interior household
budget, borrowing limit and distribution without lying about GE convergence.
"""
function stage8_independent_micro_stats(w::FMR2020.TransitionWorkspace{T},
                                        ctx::Stage3Context{T}) where {T<:AbstractFloat}
    inp=ctx.economy
    n=length(inp.B)
    gr=inp.grids.fine
    mmax=zero(T)
    weighted=zero(T)
    totalmass=zero(T)
    mincons=T(Inf)
    minborrowgap=T(Inf)
    violatingmass=zero(T)
    for t in 1:n-1,j in axes(w.distribution,2),i in eachindex(gr)
        mass=w.distribution[i,j,t]
        mass==zero(T) && continue
        gross=one(T)+w.rate[t]
        income=w.wage[t]*w.L[t]*inp.chain.levels[j]+w.profits[t]
        m=FMR2020.consumption_cost(w.terms[t],inp.params)
        consum=w.fine_consumption[i,j,t]
        saving=w.fine_savings[i,j,t]
        balance=gross*gr[i]+income-m*consum-saving
        mmax=max(mmax,abs(balance))
        weighted+=mass*abs(balance)
        totalmass+=mass
        mincons=min(mincons,consum)
        minborrowgap=min(minborrowgap,saving-inp.grids.borrowing_limit)
        if saving<inp.grids.borrowing_limit-T(1e-9) || consum<=zero(T)
            violatingmass+=mass
        end
    end
    return (max_interior_household_budget=mmax,
            mean_interior_weighted_budget=weighted/totalmass,
            minimum_occupied_consumption=mincons,
            minimum_occupied_borrowing_gap=minborrowgap,
            total_positive_mass_with_borrow_or_consumption_violation=violatingmass)
end

"""
    stage8_one_residual(ctx, candidate; partials=true)

Evaluate the exact polished strict-portfolio residual map that the Stage 5
strict solution uses. On economic infeasibility, return its stage and error.
A failed candidate is NOT replaced with clipped wealth or a solver root.
"""
function stage8_one_residual(ctx::Stage3Context{T},candidate::Vector{T};
                             partials::Bool=true) where {T<:AbstractFloat}
    n=length(ctx.economy.B)
    length(candidate)==3n-1 || throw(DimensionMismatch("599 candidate required"))
    w=FMR2020.TransitionWorkspace(ctx.economy)
    res=zeros(T,3n-1)
    values=Pair{String,Any}[]
    stage="firm prices"
    try
        FMR2020.stage3_flexible_full_prices!(w,ctx,candidate)
        append!(values,[
            "candidate_x1"=>w.terms[1],
            "candidate_K2"=>w.K[2],
            "candidate_L1"=>w.L[1],
            "candidate_first_goods_available"=>
                w.Y[1]-w.I[1]-w.adjustment[1]-
                ctx.economy.reference.export_shift*
                w.terms[1]^(-ctx.economy.params.theta_star),
        ])
    catch err
        push!(values,"failure_stage"=>stage)
        push!(values,"failure_type"=>string(typeof(err)))
        push!(values,"failure_message"=>sprint(showerror,err))
        return (ok=false,metrics=values,residual=res,work=w)
    end
    support=stage8_portfolio_support(w,ctx,w.terms[1])
    for (key,val) in pairs(support)
        push!(values,"source_guess_"*string(key)=>val)
    end
    if support.violating_mass>0 || support.offgrid_mass>0
        push!(values,"failure_stage"=>"original-state household feasibility")
        push!(values,"failure_type"=>"PositiveMassEconomicInfeasibility")
        push!(values,"failure_message"=>
            "Strict-portfolio source guess cannot fund all occupied households or leaves grid; do not run full GE from this point")
        return (ok=false,metrics=values,residual=res,work=w)
    end
    stage="full 599 economic residual"
    try
        FMR2020.stage3_flexible_full_residual!(
            res,candidate,ctx,w;
            polish_households=true,strict_portfolio=true)
    catch err
        push!(values,"failure_stage"=>stage)
        push!(values,"failure_type"=>string(typeof(err)))
        push!(values,"failure_message"=>sprint(showerror,err))
        return (ok=false,metrics=values,residual=res,work=w)
    end
    asset=view(res,1:n-1)
    labor=view(res,n:2n-1)
    goods=view(res,2n:3n-1)
    for (name,block) in (("asset",asset),("labor_with_terminal",labor),
                         ("goods_with_terminal",goods))
        idx=argmax(abs.(block))
        push!(values,"max_scaled_"*name=>maximum(abs,block))
        push!(values,"max_scaled_"*name*"_index"=>idx)
    end
    push!(values,"all599_max_scaled"=>maximum(abs,res))
    push!(values,"first_goods_scaled"=>goods[1])
    push!(values,"first_goods_level"=>goods[1]*max(one(T),ctx.economy.reference.Y))
    push!(values,"first_labor_scaled"=>labor[1])
    push!(values,"terminal_L_anchor_scaled"=>labor[end])
    push!(values,"terminal_terms_anchor_scaled"=>goods[end])
    push!(values,"first_household_C_H"=>w.aggregate_home[1])
    push!(values,"first_price_x"=>w.terms[1])
    push!(values,"first_revaluation_multiplier"=>ctx.terms_pre/w.terms[1])
    push!(values,"max_distribution_mass_error"=>maximum(w.distribution_error))
    for (key,val) in pairs(stage8_independent_micro_stats(w,ctx))
        push!(values,string(key)=>val)
    end
    push!(values,"terminal_K_gap"=>w.K[end]-ctx.economy.terminal.K)
    push!(values,"terminal_L_gap"=>w.L[end]-ctx.economy.terminal.L)
    push!(values,"terminal_terms_gap"=>w.terms[end]-ctx.economy.terminal.terms)
    if partials
        # A local 3x3 Jacobian for date-one asset/labor/goods residuals,
        # holding the OTHER 596 aggregate coordinates fixed. Never confuse
        # its step with a 599-dimensional GE root or a policy comparative static.
        columns=(1,n,2n)
        names=("K2","L1","x1")
        equation_indices=(1,n,2n)
        Jfirst=zeros(T,3,3)
        valid=true
        for j in 1:3
            col=columns[j]
            h=T(1e-5)*max(one(T),abs(candidate[col]))
            perturbed=copy(candidate)
            perturbed[col]+=h
            tmp=zeros(T,3n-1)
            ww=FMR2020.TransitionWorkspace(ctx.economy)
            try
                FMR2020.stage3_flexible_full_residual!(
                    tmp,perturbed,ctx,ww;
                    polish_households=true,strict_portfolio=true)
                for i in 1:3
                    Jfirst[i,j]=(tmp[equation_indices[i]]-
                                 res[equation_indices[i]])/h
                end
                push!(values,"dRgoods1_d"*names[j]=>Jfirst[3,j])
                push!(values,"dRasset1_d"*names[j]=>Jfirst[1,j])
                push!(values,"dRlabor1_d"*names[j]=>Jfirst[2,j])
            catch err
                valid=false
                push!(values,"partial_"*names[j]*"_status"=>
                      "UNDEFINED: "*sprint(showerror,err))
            end
        end
        if valid && all(isfinite,Jfirst)
            firstres=T[res[1],res[n],res[2n]]
            push!(values,"first_three_norm_before"=>norm(firstres))
            push!(values,"local_three_condition_number"=>cond(Jfirst))
            if cond(Jfirst)<T(1e10)
                step=-(Jfirst\firstres)
                for j in 1:3
                    push!(values,"local_Newton_delta_"*names[j]=>step[j])
                end
                # Measure ONE damped step, not an iterative/hidden solver.
                # The entire 599-vector and household distribution are
                # recomputed at the trial candidate.
                recorded=false
                for alpha in (one(T),T(0.5),T(0.25),T(0.125),T(0.0625))
                    trial=copy(candidate)
                    for j in 1:3
                        trial[columns[j]]+=alpha*step[j]
                    end
                    if !(all(isfinite,trial) &&
                         all(>(zero(T)),trial) &&
                         trial[2n]<ctx.terms_pre)
                        continue
                    end
                    rr=zeros(T,3n-1)
                    tw=FMR2020.TransitionWorkspace(ctx.economy)
                    try
                        FMR2020.stage3_flexible_full_residual!(
                            rr,trial,ctx,tw;
                            polish_households=true,strict_portfolio=true)
                    catch err
                        push!(values,"local_step_"*string(alpha)*"_infeasible"=>
                              sprint(showerror,err))
                        try
                            # The price block has ordinarily completed
                            # before an occupied-state EGM error. Distinguish
                            # genuine original wealth-budget insolvency from
                            # merely undefined interpolation.
                            ss=stage8_portfolio_support(tw,ctx,trial[2n])
                            push!(values,"local_step_"*string(alpha)*
                                  "_minimum_original_resources"=>ss.minimum_resources)
                            push!(values,"local_step_"*string(alpha)*
                                  "_violating_original_mass"=>ss.violating_mass)
                            push!(values,"local_step_"*string(alpha)*
                                  "_offgrid_original_mass"=>ss.offgrid_mass)
                        catch _
                            push!(values,"local_step_"*string(alpha)*
                                  "_feasibility_diagnostic"=>"UNAVAILABLE")
                        end
                        continue
                    end
                    newfirst=T[rr[1],rr[n],rr[2n]]
                    if norm(newfirst)<norm(firstres)
                        push!(values,"local_step_accepted"=>true)
                        push!(values,"local_step_damping"=>alpha)
                        push!(values,"local_step_x1"=>trial[2n])
                        push!(values,"local_step_K2"=>trial[1])
                        push!(values,"local_step_L1"=>trial[n])
                        push!(values,"local_first_three_norm_after"=>norm(newfirst))
                        push!(values,"local_first_asset_scaled_after"=>rr[1])
                        push!(values,"local_first_labor_scaled_after"=>rr[n])
                        push!(values,"local_first_goods_scaled_after"=>rr[2n])
                        push!(values,"local_full599_max_scaled_after"=>maximum(abs,rr))
                        push!(values,"local_first_household_C_H_after"=>tw.aggregate_home[1])
                        ss=stage8_portfolio_support(tw,ctx,trial[2n])
                        push!(values,"local_accepted_minimum_original_resources"=>
                              ss.minimum_resources)
                        push!(values,"local_accepted_violating_original_mass"=>
                              ss.violating_mass)
                        recorded=true
                        break
                    end
                end
                if !recorded
                    push!(values,"local_step_accepted"=>false)
                end
            else
                push!(values,"local_step_accepted"=>false)
                push!(values,"local_step_reason"=>"ill-conditioned three-variable local Jacobian")
            end
        else
            push!(values,"local_step_accepted"=>false)
            push!(values,"local_step_reason"=>"undefined three-variable Jacobian")
        end
    end
    return (ok=true,metrics=values,residual=res,work=w)
end

"""
    stage8_compare_strict_saved_root(ctx, path)

Read the already solved Stage 5 strict appreciation path as a diagnostic
of the SAME shock and inherited economy, never as a new numerical root.
"""
function stage8_compare_strict_saved_root(ctx::Stage3Context{T},
                                         path::AbstractString) where {T<:AbstractFloat}
    x=read_economic_root_guess(path,ctx)
    x===nothing && throw(ArgumentError("missing saved strict Stage 5 path"))
    return (impact_x=x[2*length(ctx.economy.B)],
            is_depreciation=x[2*length(ctx.economy.B)]<ctx.terms_pre)
end

"""
    stage8_read_checkpoint(path, n)

Read an accepted (but not converged) Stage 3 iterate. The checkpoint
stores ONLY x, not the Broyden Jacobian or household distribution.
The recorded max residual must be RECOMPUTED under the current strict
Euler-polished economy before a resumed full-system solve is permitted.
"""
function stage8_read_checkpoint(path::AbstractString,n::Int)
    isfile(path) || throw(ArgumentError(
        "Missing saved Stage 8 accepted-iterate checkpoint: "*path))
    x=FMR2020.stage3_read_checkpoint(path,n)
    lines=readlines(path)
    startswith(lines[3],"# iter=") ||
        throw(ArgumentError("Invalid saved iteration header"))
    startswith(lines[4],"# max_scaled_residual=") ||
        throw(ArgumentError("Invalid stored maximum-residual header"))
    prior_iter=parse(Int,split(lines[3],"=",limit=2)[2])
    prior_max=parse(Float64,split(lines[4],"=",limit=2)[2])
    prior_iter>0 && isfinite(prior_max) && prior_max>=0 ||
        throw(ArgumentError("Invalid prior checkpoint iteration or residual"))
    all(v->isfinite(v)&&v>0,x) ||
        throw(DomainError(x,"Checkpoint has a nonpositive GE unknown"))
    return (x=x,iterations=prior_iter,stored_max=prior_max)
end

"""
    stage8_verify_resume_residual(stored, actual)

Fail closed when the retained path does not satisfy the SAME full
599-dimensional economic residual map as the checkpoint source.
Checking its numerical residual is a stronger guard than trusting
a previously saved solver flag or the original source name.
"""
function stage8_verify_resume_residual(stored::Real,actual::Real)
    isfinite(stored) && isfinite(actual) && stored>=0 && actual>=0 ||
        throw(ArgumentError("Nonfinite or negative checkpoint residual"))
    isapprox(actual,stored;atol=1e-8,rtol=1e-5) ||
        throw(ErrorException(
            "Stage 8 checkpoint and CURRENT strict GE residual differ: "*
            "recorded=$(stored), re-evaluated=$(actual). "*
            "Refusing to continue a stale or different model checkpoint."))
    return true
end

"""
    stage8_append_accepted_step(path, mode, previous_iter, local_iter,
                                x, f, terms_pre, n)

Write ONLY solver-accepted iterates. The trace is a numerical log,
not a convergence or model-equilibrium certification.
"""
function stage8_append_accepted_step(path::AbstractString,
                                     mode::AbstractString,
                                     previous_iter::Int,
                                     local_iter::Int,
                                     x::AbstractVector{T},
                                     f::AbstractVector{T},
                                     terms_pre::T,
                                     n::Int) where {T<:AbstractFloat}
    length(x)==length(f)==3n-1 ||
        throw(DimensionMismatch("Invalid accepted 599-system step"))
    isfinite(terms_pre) && terms_pre>zero(T) ||
        throw(DomainError(terms_pre,"Invalid pre-shock terms"))
    all(isfinite,x) && all(isfinite,f) ||
        throw(ArgumentError("Accepted step contains nonfinite values"))
    mkpath(dirname(path))
    fresh=!isfile(path)
    open(path,"a") do io
        if fresh
            println(io,"mode,prior_accepted_iterations,local_iteration,",
                    "total_accepted_iterations,impact_x,",
                    "nominal_depreciation_logpoints,",
                    "is_depreciation,K2,L1,max_scaled_residual,",
                    "residual_norm,worst_equation_index")
        end
        t=2n
        println(io,join((
            mode,previous_iter,local_iter,
            previous_iter+local_iter,
            x[t],100*log(terms_pre/x[t]),x[t]<terms_pre,
            x[1],x[n],maximum(abs,f),norm(f),
            argmax(abs.(f))
        ),","))
    end
    return path
end

function stage8_main()
    (isempty(ARGS) || (length(ARGS)==1 && ARGS[1] in ("--solve","--resume"))) ||
        error("Usage: julia --project=. scripts/stage8_depreciation_continuation.jl [--solve|--resume]")
    mode=isempty(ARGS) ? "--audit" : ARGS[1]
    solve=mode!="--audit"
    base=normpath(joinpath(@__DIR__,".."))
    source=joinpath(base,"results","author_reference")
    roots=joinpath(base,"results","economic_reproduction")
    output=joinpath(base,"results","stage8_depreciation")
    checkpoint=joinpath(output,"source_initialized_UNCONVERGED.checkpoint")
    trace=joinpath(output,"accepted_full599_trace.csv")
    steady=FMR2020.solve_three_steady_states(FMR2020.default_parameters();
                                             verbose=false)
    boom=stage5_existing_boom(steady,roots)
    ctx=FMR2020.stage3_reversal_context(
        steady,boom;shock_start=41,leverage_hat=1.0/16.0)
    n=length(ctx.economy.B)
    if mode=="--resume"
        # Do NOT repeat the 399-source audit or the initial 3x3 trial.
        # Restart from the ACTUAL accepted full-599 iterate on the PC.
        saved=stage8_read_checkpoint(checkpoint,3n-1)
        trial=stage8_one_residual(ctx,saved.x;partials=false)
        trial.ok || error(
            "Saved accepted-iterate checkpoint is infeasible/undefined under "*
            "CURRENT strict 599 equations. No resume. "*string(trial.metrics))
        now=maximum(abs,trial.residual)
        stage8_verify_resume_residual(saved.stored_max,now)
        values=Pair{String,Any}[
            "model"=>"strict-portfolio, Euler-polished Julia boom initial state",
            "resume_from_accepted_iteration"=>saved.iterations,
            "checkpoint_max_scaled_recorded"=>saved.stored_max,
            "checkpoint_max_scaled_recomputed"=>now,
            "checkpoint_impact_x"=>saved.x[2n],
            "checkpoint_nominal_depreciation_logpoints"=>
                100*log(ctx.terms_pre/saved.x[2n]),
            "checkpoint_is_depreciation"=>saved.x[2n]<ctx.terms_pre,
            "equilibrium_root_status"=>"NOT CERTIFIED",
        ]
        append!(values,trial.metrics)
        stage8_store_metrics(joinpath(output,"resume_candidate_residuals.csv"),
                             values)
        println("STAGE8 VERIFIED CHECKPOINT: accepted previous iterations=",
                saved.iterations," current strict max=",now,
                " x1=",saved.x[2n],
                " depreciation=",saved.x[2n]<ctx.terms_pre)
        if saved.x[2n]>=ctx.terms_pre
            println("CAUTION: accepted iterate has crossed into appreciation. ",
                    "Its eventual convergence would not establish a second ",
                    "depreciation equilibrium.")
        end
        candidate=saved.x
        previous_iter=saved.iterations
    else
        candidate,meta=stage8_author_initial_guess(source,ctx)
        values=Pair{String,Any}[
            "model"=>"strict-portfolio, Euler-polished Julia boom initial state",
            "source_guess_only"=>true,
            "equilibrium_root_status"=>"NOT CERTIFIED",
            "unknown_count"=>length(candidate),
            "asset_market_equations"=>n-1,
            "labor_market_equations"=>n-1,
            "goods_market_equations"=>n-1,
            "terminal_anchors"=>2,
            "original_source_K_initial"=>meta[1],
            "current_strict_K_initial"=>ctx.economy.initial.K,
            "original_source_q_initial"=>meta[2],
            "current_strict_q_initial"=>ctx.economy.q_initial,
            "original_source_B_initial"=>meta[3],
            "current_strict_B_initial"=>ctx.economy.B_initial,
            "original_source_terms_pre"=>meta[4],
            "current_strict_terms_pre"=>ctx.terms_pre,
            "foreign_credit_B1"=>ctx.economy.B[1],
            "foreign_credit_Bterminal"=>ctx.economy.B[end],
        ]
        strict_saved=joinpath(
            roots,"economic_flexible_leveraged_euler_polished_"*
            "STRICT_PORTFOLIO_MARKET_ROOT_UNCERTIFIED.csv")
        if isfile(strict_saved)
            high=stage8_compare_strict_saved_root(ctx,strict_saved)
            push!(values,"saved_strict_stage5_root_available"=>true)
            push!(values,"prior_strict_root_impact_x"=>high.impact_x)
            push!(values,"prior_strict_root_depreciates"=>high.is_depreciation)
        else
            push!(values,"saved_strict_stage5_root_available"=>false)
            push!(values,"prior_strict_root_impact_x"=>"NOT VERIFIED IN THIS RUN")
        end
        baseline=stage8_one_residual(ctx,candidate)
        append!(values,baseline.metrics)
        push!(values,"source_candidate_599_status"=>
              (baseline.ok ? "DEFINED_NOT_CERTIFIED" :
                             "UNDEFINED_ECONOMIC_CANDIDATE"))
        stage8_store_metrics(joinpath(output,"initial_candidate_residuals.csv"),
                             values)
        for (k,v) in values
            println(k," = ",v)
        end
        if !solve
            println("STAGE8 AUDIT ONLY: no 599-variable GE solve")
            return nothing
        end
        baseline.ok || error("Source depreciation guess economically undefined")
        # Do not silently destroy an existing accepted-iterate checkpoint.
        !isfile(checkpoint) || error(
            "Stage 8 already has an accepted-iterate checkpoint; "*
            "use --resume rather than --solve to preserve progress.")
        previous_iter=0
    end

    # Numerical restart, NOT a warm-start of the Broyden matrix.
    # The existing solver rebuilds its full Jacobian at resumed iteration 0.
    # No implementation change or novel 599-dimensional solver is needed.
    last=(x,f,i)->begin
        stage8_append_accepted_step(
            trace,mode,previous_iter,i,x,f,ctx.terms_pre,n)
        FMR2020.stage3_write_checkpoint(
            checkpoint,x,f,previous_iter+i)
        println("STAGE8 ACCEPTED total=",previous_iter+i,
                " x1=",x[2n],
                " depreciation=",x[2n]<ctx.terms_pre,
                " max_scaled=",maximum(abs,f))
    end
    maxiter=mode=="--resume" ? 12 : 5
    maxeval=mode=="--resume" ? 2200 : 1000
    println("STARTING STRICT 599 FULL GE mode=",mode,
            " from accepted=",previous_iter,
            " maxiter=",maxiter," max_evaluations=",maxeval)
    sol=FMR2020.solve_stage3_flexible_full(
        ctx;x0=candidate,polish_households=true,strict_portfolio=true,
        tolerance=1e-6,maxiter=maxiter,max_evaluations=maxeval,
        on_accepted_step=last,verbose=true)
    sol.converged || error("Stage 8 did not converge: no equilibrium claim")
    diagnostics=FMR2020.stage3_diagnostics(
        sol,ctx;explicit_flexible=true)
    certified_market=(
        diagnostics.maximum_scaled_solver_residual<1e-6 &&
        diagnostics.max_household_budget<1e-8 &&
        diagnostics.max_mass_deviation<1e-8 &&
        diagnostics.max_asset_market_level<1e-4 &&
        diagnostics.max_labor_market_level<1e-4 &&
        diagnostics.max_goods_market_level<1e-4)
    certified_market || error(
        "Returned Stage 8 root fails independent market/budget diagnostics")
    marketpath=joinpath(
        output,"source_initialized_STRICT_MARKET_ROOT_UNCERTIFIED.csv")
    FMR2020.stage3_write_csv(sol,ctx,marketpath)
    println("STAGE8 SAVED MARKET ROOT (NOT CERTIFIED INFINITE HORIZON): ",
            marketpath)
    println("SOURCE-INITIALIZED ROOT PRICE X1=",sol.workspace.terms[1],
            " prior Xpre=",ctx.terms_pre,
            " nominal_depreciation_logpoints=",
            100*log(ctx.terms_pre/sol.workspace.terms[1]))
    println("INDEPENDENT DIAGNOSTICS ",diagnostics)
    return nothing
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    stage8_main()
end
