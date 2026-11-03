# Evaluate the author's ACTUAL MATLAB 399-dimensional saved paths
# against the explicit flexible-price 599-dimensional economic closure.
# Unlike evaluate_author_matlab_solution.jl, this script does not re-solve
# the boom or the nested 399 system; all shock-date inheritance data are
# read from the source .mat-exported fixture.
#
# Run from FMR-2020:
# julia --project=. scripts/evaluate_author_explicit_goods.jl results/author_reference
# Optional expensive full nonlinear solve, initialized at the author path:
# julia --project=. scripts/evaluate_author_explicit_goods.jl results/author_reference --solve

using FMR2020
using LinearAlgebra

function source_vector(dir::AbstractString, label::AbstractString)
    file=joinpath(dir,string(label,".csv"))
    isfile(file) || error("Missing actual source MATLAB CSV: $file")
    values=parse.(Float64,strip.(readlines(file)))
    all(isfinite,values) || error("Nonfinite actual source MATLAB CSV: $file")
    return values
end

function source_matrix(dir::AbstractString, label::AbstractString)
    file=joinpath(dir,string(label,".csv"))
    isfile(file) || error("Missing actual source MATLAB CSV: $file")
    data=readlines(file)
    isempty(data) && error("Empty actual source MATLAB CSV: $file")
    ncol=length(split(strip(first(data)),","))
    result=Matrix{Float64}(undef,length(data),ncol)
    for (i,line) in enumerate(data)
        cols=split(strip(line),",")
        length(cols)==ncol || error("Ragged actual source MATLAB CSV: $file")
        for j in 1:ncol
            result[i,j]=parse(Float64,strip(cols[j]))
        end
    end
    all(isfinite,result) || error("Nonfinite actual source MATLAB CSV: $file")
    return result
end

function source_reversal_context(steady, mu::Matrix{Float64},
                                 metadata::Vector{Float64},
                                 leverage::Float64)
    length(metadata)==7 || error("Expected seven original MATLAB inherited scalars")
    n=steady.parameters.transition_horizon
    n==200 || error("Author sudden-stop T=200 required")
    gr=FMR2020.stage3_full_grids(steady.grids)
    size(mu)==(length(gr.fine),length(steady.chain.levels)) ||
        throw(DimensionMismatch("Actual MATLAB shock-date distribution has wrong shape"))
    abs(sum(mu)-1)<1e-8 ||
        error("Actual author MATLAB shock-date distribution mass is not one")
    abs(steady.calibrated.K-metadata[5])<1e-4 ||
        error("Julia calibrated terminal capital differs materially from original MATLAB")
    abs(steady.calibrated.beta-metadata[6])<1e-6 ||
        error("Julia calibrated beta differs materially from original MATLAB")
    initial=FMR2020.stage3_initial_state(
        steady.initial,metadata[1],steady.initial.L,metadata[3],
        metadata[2],metadata[4],mu)
    # Exactly the author's first-quarter B shock and quartic reversal:
    # Bnew(t)=Bcal+(Bprev-Bcal)*(1-(t-1)/(T-1))^4.
    Bcal=steady.calibrated.B
    Bnew=[Bcal+(metadata[3]-Bcal)*(1-(t-1)/(n-1))^4
          for t in 1:n]
    economy=FMR2020.TransitionInputs(
        steady.parameters,steady.reference,steady.chain,gr,
        initial,steady.calibrated,Bnew,ones(Float64,n),17.0,
        metadata[2],metadata[3])
    return FMR2020.Stage3Context(economy,copy(mu),metadata[4],
                                 leverage,100.0,true,
                                 "author saved path, explicit full goods")
end

function inspect_author_path(steady,dir::String,mu::Matrix{Float64},
                             metadata::Vector{Float64},
                             name::String,lev::Float64;
                             solve::Bool=false)
    ctx=source_reversal_context(steady,mu,metadata,lev)
    original=source_vector(dir,name*"_candidate")
    original_terms=source_vector(dir,name*"_source_terms")
    original_outer=source_vector(dir,name*"_source_residual")
    n=length(ctx.economy.B)
    length(original)==2n-1 && length(original_outer)==2n-1 &&
        length(original_terms)==n ||
        throw(DimensionMismatch("Author original saved 399/200 path dimensions"))
    x=vcat(original,original_terms)
    w=FMR2020.TransitionWorkspace(ctx.economy)
    res=zeros(Float64,3n-1)
    println("AUTHOR SAVED ",name,
            " MATLAB raw 399 max=",maximum(abs,original_outer),
            " (NOT evidence of full goods clearing)")
    try
        FMR2020.stage3_flexible_full_residual!(res,x,ctx,w)
    catch err
        println("AUTHOR SAVED ",name,
                " full Julia 599 residual UNDEFINED: ",
                sprint(showerror,err))
        println("Source 399 root cannot be certified as strict-budget Julia GE.")
        return false
    end
    coremax=maximum(abs,view(res,1:2n-1))
    goods=view(res,2n:3n-1)
    maxgoods=maximum(abs,goods)
    firstgoods=goods[1]
    println("AUTHOR SAVED ",name,
            " explicit 599 evaluation: max_KL_scaled=",coremax,
            " first_goods_scaled=",firstgoods,
            " max_goods_scaled=",maxgoods,
            " all_599_max_scaled=",maximum(abs,res))
    if !solve
        println("EVALUATION ONLY. A nonzero 599-vector is NOT a GE.")
        return coremax<1e-6 && maxgoods<1e-6
    end
    println("STARTING FULL 599 NONLINEAR SOLVE from ACTUAL author MATLAB path")
    sol=FMR2020.solve_stage3_flexible_full(
        ctx;x0=x,tolerance=1e-6,maxiter=60,
        max_evaluations=5000,verbose=true)
    sol.converged || error("Explicit 599 full GE returned unconverged")
    diag=FMR2020.stage3_diagnostics(sol,ctx;fixed=false,
                                   explicit_flexible=true)
    println("EXPLICIT FLEX GE CANDIDATE ",name," ",diag)
    (diag.max_asset_market_level<1e-4 &&
     diag.max_labor_market_level<1e-4 &&
     diag.max_goods_market_level<1e-4 &&
     diag.max_household_budget<1e-8 &&
     diag.max_mass_deviation<1e-8 &&
     diag.max_dynamic_euler<1e-4 &&
     diag.weighted_mean_dynamic_euler<1e-5) ||
        error("Explicit 599 candidate failed independent economic checks")
    output=joinpath(dir,"julia_full_599_"*name*".csv")
    FMR2020.stage3_write_csv(sol,ctx,output)
    println("FULL 599 GE VERIFIED AND SAVED: ",output,
            " (including Euler/KKT tolerances)")
    return true
end

function main()
    (length(ARGS)==1 ||
     (length(ARGS)==2 &&
      ARGS[2] in ("--solve", "--solve-unleveraged", "--solve-leveraged"))) ||
        error("Usage: julia --project=. scripts/evaluate_author_explicit_goods.jl FIXTURE_DIR [--solve|--solve-unleveraged|--solve-leveraged]")
    dir=ARGS[1]
    mode=length(ARGS)==1 ? "--audit" : ARGS[2]
    do_solve=mode!="--audit"
    st=FMR2020.solve_three_steady_states(FMR2020.default_parameters();
                                           verbose=false)
    mu=source_matrix(dir,"shock_pre_distribution")
    meta=source_vector(dir,"shock_inherited_metadata")
    states=mode=="--solve-unleveraged" ? [("unleveraged",0.0)] :
           mode=="--solve-leveraged" ? [("leveraged",1.0/16.0)] :
           [("unleveraged",0.0),("leveraged",1.0/16.0)]
    results=Bool[]
    for (name,lev) in states
        push!(results,inspect_author_path(st,dir,mu,meta,name,lev;
                                          solve=do_solve))
    end
    if do_solve
        all(results) ||
            error("Explicit 599 GE not verified. Do not claim convergence.")
    else
        println("SOURCE-PATH AUDIT COMPLETE: ",count(identity,results),
                "/",length(results)," original saved paths satisfy the " *
                "strict explicit 599 GE. Nonzero goods residuals are " *
                "a measured finding, not a file or program error.")
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end
