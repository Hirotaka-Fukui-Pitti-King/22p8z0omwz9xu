# Export the actual shock-date distribution and both trade-price
# branches. MATLAB/Octave revalues EXACTLY the same positive-mass
# household portfolio nodes using original HANKSOME basefun.m.
# No 399-variable sudden-stop nonlinear solve is performed.
using FMR2020

function savecsv(path::AbstractString,x::AbstractArray)
    mat=ndims(x)==1 ? reshape(x,:,1) : x
    ndims(mat)==2 || throw(DimensionMismatch("expected vector/matrix"))
    open(path,"w") do io
        for i in axes(mat,1)
            println(io,join((repr(mat[i,j]) for j in axes(mat,2)),","))
        end
    end
end

function main()
    (length(ARGS)==1 ||
     (length(ARGS)==2 && ARGS[2]=="--checkpoint")) ||
        error("Usage: julia export_stage3_fx_fixture.jl OUTPUT_DIR [--checkpoint]")
    directory=ARGS[1]
    checkpoint_mode=length(ARGS)==2
    mkpath(directory)
    root=normpath(joinpath(@__DIR__,".."))
    st=solve_three_steady_states(default_parameters();verbose=false)
    boomctx=stage3_expansion_context(st)
    csvfile=joinpath(root,"results","stage2_transition_expansion.csv")
    isfile(csvfile) || error("Stage2 GE starting-path CSV required")
    lines=readlines(csvfile)
    n=length(boomctx.economy.B)
    length(lines)==n+1 || error("Stage2 CSV transition horizon mismatch")
    names=split(lines[1],",")
    ki=findfirst(==("K"),names)
    li=findfirst(==("L"),names)
    (isnothing(ki) || isnothing(li)) &&
        error("Stage2 CSV missing K/L columns")
    k=Float64[]
    l=Float64[]
    for line in lines[2:end]
        cells=split(line,",")
        push!(k,parse(Float64,cells[ki]))
        push!(l,parse(Float64,cells[li]))
    end
    guess=vcat(k[2:end],l)
    boom=solve_stage3_flexible(boomctx;x0=guess,
                               tolerance=1e-6,maxiter=3,verbose=false)
    boom.converged || error("Pre-shock foreign-denominated boom GE is not solved")
    ctx=stage3_reversal_context(st,boom;shock_start=41,
                                leverage_hat=0.0)
    inp=ctx.economy
    x=zeros(Float64,2*n-1)
    for t in 1:n-1
        x[t]=st.calibrated.K+(inp.initial.K-st.calibrated.K)*0.8^t
    end
    x[n:end].=st.calibrated.L
    # Isolate the last ACCEPTED but UNCONVERGED 399-variable iterate.
    # Never mistake it for a solved GE. This is an especially important
    # off-equilibrium comparison: the original source-like EGM map
    # may deviate from production much more here than at the first guess.
    if checkpoint_mode
        saved=joinpath(root,"results",
                       "stage3_unleveraged_UNCONVERGED.checkpoint")
        isfile(saved) ||
            error("No unleveraged accepted-iterate checkpoint at $saved")
        x=FMR2020.stage3_read_checkpoint(saved,2*n-1)
        println("STAGE3 SOURCE-MAP AUDIT: accepted but UNCONVERGED " *
                "checkpoint; 399 residuals will be reevaluated.")
    else
        println("STAGE3 SOURCE-MAP AUDIT: original initial K/L guess.")
    end
    # Export ALL inputs needed to call the author's unmodified
    # Transition/solve_transition_bis.m at the exact same 399-vector.
    # This is a direct original-source residual evaluation, NOT fsolve.
    p=inp.params
    savecsv(joinpath(directory,"candidate_399.csv"),x)
    savecsv(joinpath(directory,"B_stop.csv"),inp.B)
    savecsv(joinpath(directory,"z_path.csv"),inp.z)
    savecsv(joinpath(directory,"coarse_grid.csv"),st.grids.coarse_full)
    savecsv(joinpath(directory,"productivity.csv"),inp.chain.levels)
    savecsv(joinpath(directory,"transition.csv"),inp.chain.transition)
    savecsv(joinpath(directory,"source_parameters.csv"),Float64[
        st.calibrated.beta,p.chi,p.theta,p.gamma,p.alpha,p.delta,
        p.epsilon_f,p.epsilon_w,inp.phi_ac,p.theta_star,
        inp.reference.export_shift,inp.reference.chi_dis,
        inp.terminal.K,p.grid_min_z,p.grid_power])
    grid=inp.grids.fine
    savecsv(joinpath(directory,"fine_grid.csv"),grid)
    savecsv(joinpath(directory,"mu_pre.csv"),ctx.mu_pre)
    meta=zeros(Float64,2,11)
    for (row,lev) in enumerate((0.0,1/16))
        case=stage3_reversal_context(st,boom;shock_start=41,
                                     leverage_hat=lev)
        w=FMR2020.TransitionWorkspace(case.economy)
        FMR2020.stage3_flexible_prices!(w,case,x)
        p=case.economy.params
        budget_out=case.economy.B[1]+w.q[1]*w.K[2]
        other=w.wage[1]*w.L[1]+w.profits[1]
        gross=1+w.rate[1]
        C=w.Y[1]-w.I[1]-w.adjustment[1]-
          p.chi*(gross*case.economy.q_initial*w.K[1]-budget_out+other)
        E=-p.chi*gross*case.economy.B_initial*case.terms_pre
        D=case.economy.reference.export_shift
        roots=FMR2020.foreign_debt_trade_roots(C,E,D,p.theta_star)
        roots.high===nothing && error("Expected double-root shock candidate")
        meta[row,:].=[w.K[1],case.economy.q_initial,
                      case.economy.B_initial,case.terms_pre,
                      lev,C,E,D,p.theta_star,roots.low,roots.high]
        # The unleveraged case selects the SAME low trade root in
        # original source and current Julia. Hold the entire candidate
        # K/L and price path fixed and compare all 199 backward EGM steps.
        # This identifies a policy-map difference without root-finding.
        if lev==0.0
            FMR2020.stage3_household_backward!(w,case)
            # Evaluate occupied probabilities under the ORIGINAL
            # Julia candidate path, without claiming it is a GE.
            FMR2020.stage3_distribution_forward!(w,case)
            # Store the unscaled Julia 399-vector for DIRECT comparison
            # with the original MATLAB source's raw residual vector.
            prod=zeros(Float64,2*n-1)
            FMR2020.stage3_market_residual!(prod,w,case)
            # Numerical identity check for the NEW 599-flex residual map:
            # supplying the existing 399-map terms of trade as EXPLICIT
            # unknowns must reproduce ALL 399 asset/labor residuals.
            # Unlike the nested original, the last 200 new residuals
            # are ACTUAL household goods-market clearing, not an
            # assumed scalar fsolve root. This is a map test, NOT GE.
            explicit_x=vcat(x,copy(w.terms))
            explicit_res=zeros(Float64,3*n-1)
            explicit_work=FMR2020.TransitionWorkspace(case.economy)
            FMR2020.stage3_flexible_full_residual!(
                explicit_res,explicit_x,case,explicit_work)
            error399=maximum(abs.(explicit_res[1:2*n-1].-prod))
            error399<1e-8 ||
                error("Explicit flexible map's 399 K/L conditions differ: $error399")
            @assert all(isfinite,explicit_res)
            goods_max=maximum(abs,view(explicit_res,2*n:3*n-1))
            println("EXPLICIT 599 FLEX NUMERICAL MAP PASS: " *
                    "legacy_KL_maxdifference=",error399,
                    " actual_household_goods_maxscaled=",goods_max,
                    " new_residual_dimension=",length(explicit_res),
                    "; candidate is NOT a solved GE.")
            raw=similar(prod)
            stage3_unscale_residual!(raw,prod,case)
            savecsv(joinpath(directory,"julia_raw_399.csv"),raw)
            savecsv(joinpath(directory,"julia_terms.csv"),w.terms)
            savecsv(joinpath(directory,"julia_rate.csv"),w.rate)
            savecsv(joinpath(directory,"julia_wage.csv"),w.wage)
            savecsv(joinpath(directory,"julia_aggregate_C.csv"),
                    vcat(w.aggregate_home[1:end-1],w.aggregate_home[end-1]))
            savecsv(joinpath(directory,"julia_aggregate_A.csv"),
                    w.aggregate_savings[1:end-1])
            calibrated=st.calibrated
            refpolicy=source_matlab_stationary_egm(
                st.grids,st.chain,p,calibrated.beta,calibrated.r,
                calibrated.terms,calibrated.wage,calibrated.L,
                calibrated.profits;tolerance=1e-9,maxiter=5000)
            future=refpolicy.coarse_consumption
            nf=length(case.economy.grids.fine)
            ns=length(case.economy.chain.levels)
            source_fine_c=Array{Float64}(undef,nf,ns,n)
            source_fine_a=Array{Float64}(undef,nf,ns,n)
            source_fine_c[:,:,n].=refpolicy.fine_consumption
            source_fine_a[:,:,n].=refpolicy.fine_savings
            comparison=joinpath(directory,"unleveraged_backward_egm.csv")
            open(comparison,"w") do io
                println(io,"t,max_coarse_abs_diff,max_fine_defined_abs_diff," *
                           "reference_fine_nan,production_fine_nan," *
                           "occupied_source_nan_mass,max_occupied_source_budget_gap," *
                           "max_occupied_fine_abs_diff,weighted_occupied_fine_abs_diff")
                for t in (n-1):-1:1
                    step=source_matlab_egm_step(
                        st.grids,st.chain,p,future,calibrated.beta,
                        w.rate[t],w.rate[t+1],w.terms[t],w.wage[t],
                        w.L[t],w.profits[t])
                    source_fine_c[:,:,t].=step.fine_consumption
                    source_fine_a[:,:,t].=step.fine_savings
                    oldcoarse=view(w.coarse_consumption,:,:,t)
                    sourceactive=view(step.coarse_consumption,33:193,:)
                    maxcoarse=maximum(abs.(sourceactive.-oldcoarse))
                    oldfine=view(w.fine_consumption,:,:,t)
                    mask=isfinite.(step.fine_consumption) .&
                         isfinite.(oldfine)
                    any(mask) || error("No jointly defined fine policies at time $t")
                    maxfine=maximum(abs.(step.fine_consumption[mask].-
                                         oldfine[mask]))
                    invalidmass=0.0
                    maxbudget=0.0
                    maxoccupied=0.0
                    meanoccupied=0.0
                    distribution=view(w.distribution,:,:,t)
                    for j in axes(distribution,2),i in axes(distribution,1)
                        mass=distribution[i,j]
                        mass>0 || continue
                        if !isfinite(step.fine_consumption[i,j]) ||
                           !isfinite(step.fine_savings[i,j])
                            invalidmass+=mass
                        else
                            dc=abs(step.fine_consumption[i,j]-oldfine[i,j])
                            maxoccupied=max(maxoccupied,dc)
                            meanoccupied+=mass*dc
                            maxbudget=max(maxbudget,abs(step.budget_gap[i,j]))
                        end
                    end
                    println(io,join((t,maxcoarse,maxfine,
                                     count(isnan,step.fine_consumption),
                                     count(isnan,oldfine),
                                     invalidmass,maxbudget,
                                     maxoccupied,meanoccupied),","))
                    if t in (1,2,10,50,100,n-1)
                        println("SOURCE vs PRODUCTION BACKWARD EGM t=",t,
                                " coarse_max=",maxcoarse,
                                " fine_defined_max=",maxfine,
                                " source_nan=",count(isnan,step.fine_consumption),
                                " occupied_source_nan_mass=",invalidmass,
                                " occupied_budget_gap=",maxbudget,
                                " occupied_max_fine_diff=",maxoccupied,
                                " occupied_weighted_fine_diff=",meanoccupied)
                    end
                    future=step.coarse_consumption
                end
            end
            println("BACKWARD EGM COMPARISON COMPLETE: ",comparison)
            # Source-style full 399-element market residual comparison.
            # Prices/K/L/initial FX mass are held fixed, while source
            # 193-node EGM and its induced 769-node distribution evolve
            # independently. The result is NOT a solved equilibrium.
            mu_source=copy(view(w.distribution,:,:,1))
            next_source=similar(mu_source)
            source_ch=zeros(Float64,n)
            source_a=zeros(Float64,n)
            first_invalid=0
            for t in 1:n
                for j in axes(mu_source,2),i in axes(mu_source,1)
                    mass=mu_source[i,j]
                    mass>0 || continue
                    cc=source_fine_c[i,j,t]
                    aa=source_fine_a[i,j,t]
                    if !(isfinite(cc) && cc>0 && isfinite(aa) &&
                         aa>=case.economy.grids.borrowing_limit)
                        first_invalid=t
                        break
                    end
                    source_ch[t]+=mass*cc
                    source_a[t]+=mass*aa
                end
                first_invalid>0 && break
                if t<n
                    FMR2020.distribution_step!(
                        next_source,mu_source,view(source_fine_a,:,:,t),
                        case.economy.grids.fine,case.economy.chain.transition)
                    mu_source,next_source=next_source,mu_source
                end
            end
            if first_invalid>0
                println("SOURCE 399-MAP INFEASIBLE: positive-mass original " *
                        "policy at t=",first_invalid,
                        "; cannot claim an economically valid GE")
            else
                scaleK=max(1.0,case.economy.initial.K)
                scaleL=max(1.0,case.economy.initial.L)
                unionfactor=p.chi*(p.epsilon_w-1)/p.epsilon_w
                ra=zeros(Float64,n-1)
                rl=zeros(Float64,n)
                for t in 1:n-1
                    ra[t]=(source_a[t]-w.aggregate_savings[t])/
                           w.q[t]/scaleK
                end
                # Source solve_transition_bis.m:163 extends the
                # penultimate aggregate consumption to final t=n:
                # c_t=[c_t';c_t(end)]. This matters for labor residual n.
                for t in 1:n
                    ct=t<n ? t : n-1
                    source_l=unionfactor*w.wage[t]/
                             (case.economy.reference.chi_dis*source_ch[ct])
                    production_l=unionfactor*w.wage[t]/
                             (case.economy.reference.chi_dis*w.aggregate_home[ct])
                    rl[t]=(source_l-production_l)/scaleL
                end
                maxa=maximum(abs,ra)
                maxl=maximum(abs,rl)
                ai=argmax(abs.(ra))
                li=argmax(abs.(rl))
                println("SOURCE 399 RESIDUAL MAP vs PRODUCTION SAME X " *
                        "max_asset_scaled_diff=",maxa,
                        " at_t=",ai," max_labor_scaled_diff=",maxl,
                        " at_t=",li)
                println("SOURCE 399 RESIDUAL MAP t1 source_A=",source_a[1],
                        " current_A=",w.aggregate_savings[1],
                        " source_CH=",source_ch[1],
                        " current_CH=",w.aggregate_home[1])
                open(joinpath(directory,"unleveraged_source_399_map.csv"),"w") do io
                    println(io,"t,source_CH,production_CH,source_assets,production_assets," *
                               "source_minus_production_asset_residual_scaled," *
                               "source_minus_production_labor_residual_scaled")
                    for t in 1:n
                        av=t<n ? ra[t] : NaN
                        println(io,join((t,source_ch[t],w.aggregate_home[t],
                                         source_a[t],w.aggregate_savings[t],
                                         av,rl[t]),","))
                    end
                end
            end
        end
        println("FX REAL INITIAL CASE leverage=",lev,
                " roots=(",roots.low,",",roots.high,")",
                " Julia selected=",w.terms[1])
    end
    savecsv(joinpath(directory,"metadata.csv"),meta)
    println("SHOCK FX FIXTURE INPUT EXPORTED (no sudden-stop GE solve)")
end

main()
