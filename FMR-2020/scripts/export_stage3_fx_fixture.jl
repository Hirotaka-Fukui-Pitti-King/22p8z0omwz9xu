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
    length(ARGS)==1 || error("Pass FX fixture output directory")
    directory=ARGS[1]
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
            calibrated=st.calibrated
            refpolicy=source_matlab_stationary_egm(
                st.grids,st.chain,p,calibrated.beta,calibrated.r,
                calibrated.terms,calibrated.wage,calibrated.L,
                calibrated.profits;tolerance=1e-9,maxiter=5000)
            future=refpolicy.coarse_consumption
            comparison=joinpath(directory,"unleveraged_backward_egm.csv")
            open(comparison,"w") do io
                println(io,"t,max_coarse_abs_diff,max_fine_defined_abs_diff," *
                           "reference_fine_nan,production_fine_nan," *
                           "occupied_source_nan_mass,max_occupied_source_budget_gap")
                for t in (n-1):-1:1
                    step=source_matlab_egm_step(
                        st.grids,st.chain,p,future,calibrated.beta,
                        w.rate[t],w.rate[t+1],w.terms[t],w.wage[t],
                        w.L[t],w.profits[t])
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
                    if t==1
                        initialmass=view(w.distribution,:,:,1)
                        for j in axes(initialmass,2),i in axes(initialmass,1)
                            mass=initialmass[i,j]
                            mass>0 || continue
                            if !isfinite(step.fine_consumption[i,j]) ||
                               !isfinite(step.fine_savings[i,j])
                                invalidmass+=mass
                            else
                                maxbudget=max(maxbudget,
                                              abs(step.budget_gap[i,j]))
                            end
                        end
                    end
                    println(io,join((t,maxcoarse,maxfine,
                                     count(isnan,step.fine_consumption),
                                     count(isnan,oldfine),
                                     invalidmass,maxbudget),","))
                    if t in (1,2,10,50,100,n-1)
                        println("SOURCE vs PRODUCTION BACKWARD EGM t=",t,
                                " coarse_max=",maxcoarse,
                                " fine_defined_max=",maxfine,
                                " source_nan=",count(isnan,step.fine_consumption),
                                " occupied_source_nan_mass=",invalidmass,
                                " occupied_budget_gap=",maxbudget)
                    end
                    future=step.coarse_consumption
                end
            end
            println("BACKWARD EGM COMPARISON COMPLETE: ",comparison)
        end
        println("FX REAL INITIAL CASE leverage=",lev,
                " roots=(",roots.low,",",roots.high,")",
                " Julia selected=",w.terms[1])
    end
    savecsv(joinpath(directory,"metadata.csv"),meta)
    println("SHOCK FX FIXTURE INPUT EXPORTED (no sudden-stop GE solve)")
end

main()
