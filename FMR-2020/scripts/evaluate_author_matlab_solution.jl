# Evaluate ACTUAL user-supplied author's MATLAB solved 399/599 paths.
# This is NOT a nonlinear solve. Original reference source: exp_flex.mat,
# exp_flex_lev16.mat, exp_both_41_lev16.mat (17-digit CSV fixture bundle).
# Invoke from FMR-2020 with:
# julia --project=. scripts/evaluate_author_matlab_solution.jl results/author_reference
using FMR2020
using LinearAlgebra

function reference_vector(dir::AbstractString,basename::AbstractString)
    path=joinpath(dir,string(basename,".csv"))
    isfile(path) || error("Missing MATLAB reference: $path")
    rows=readlines(path)
    x=parse.(Float64,strip.(rows))
    all(isfinite,x) || error("Nonfinite MATLAB reference: $path")
    return x
end

function reference_matrix(dir::AbstractString,basename::AbstractString)
    path=joinpath(dir,string(basename,".csv"))
    isfile(path) || error("Missing MATLAB matrix: $path")
    lines=readlines(path)
    isempty(lines) && error("Empty MATLAB matrix: $path")
    cols=split(strip(lines[1]),",")
    ncol=length(cols)
    result=Matrix{Float64}(undef,length(lines),ncol)
    for (i,line) in enumerate(lines)
        cells=split(strip(line),",")
        length(cells)==ncol || error("Inconsistent MATLAB CSV: $path")
        for j in 1:ncol
            result[i,j]=parse(Float64,strip(cells[j]))
        end
    end
    all(isfinite,result) || error("Nonfinite MATLAB matrix: $path")
    return result
end

function check_author_solved_path(dir::AbstractString,label::String,
                                  ctx::Stage3Context{Float64};
                                  fixed::Bool=false)
    x=reference_vector(dir,label*"_candidate")
    source=reference_vector(dir,label*"_source_residual")
    source_terms=reference_vector(dir,label*"_source_terms")
    n=length(ctx.economy.B)
    dim=fixed ? 3*n-1 : 2*n-1
    length(x)==dim && length(source)==dim && length(source_terms)==n ||
        throw(DimensionMismatch("Author reference path has incorrect T=200 dimensions"))
    source_max=maximum(abs,source)
    source_max<1e-6 || error("Author file is not a solved source GE: $label")
    workspace=FMR2020.TransitionWorkspace(ctx.economy)
    normalized=zeros(Float64,dim)
    println("MATLAB SOLVED ",label," source_max_raw=",source_max,
            " source_L2_raw=",norm(source))
    try
        if fixed
            FMR2020.stage3_fixed_residual!(normalized,x,ctx,workspace)
        else
            FMR2020.stage3_residual!(normalized,x,ctx,workspace)
        end
    catch err
        println("JULIA MAP UNDEFINED at ACTUAL author-solved x for ",
                label,": ",sprint(showerror,err))
        return false
    end
    unscaled=similar(normalized)
    stage3_unscale_residual!(unscaled,normalized,ctx;fixed=fixed)
    discrep=abs.(unscaled.-source)
    idx=argmax(discrep)
    terms_error=maximum(abs.(workspace.terms.-source_terms))
    println("JULIA EVALUATED ORIGINAL SOLUTION ",label,
            " scaled_max=",maximum(abs,normalized),
            " raw_max=",maximum(abs,unscaled),
            " MATLAB_vs_Julia_raw_resid_max=",discrep[idx],
            " at_raw_index=",idx,
            " term_initial_author=",source_terms[1],
            " term_initial_julia=",workspace.terms[1],
            " max_terms_error=",terms_error)
    println("FIRST FIVE RAW SOURCE RESIDUALS: ",source[1:5])
    println("FIRST FIVE RAW JULIA RESIDUALS: ",unscaled[1:5])
    return maximum(abs,normalized)<1e-6 &&
           maximum(discrep)<1e-5 &&
           terms_error<1e-5
end

function main()
    length(ARGS)==1 ||
        error("Usage: julia --project=. scripts/evaluate_author_matlab_solution.jl FIXTURE_DIR")
    dir=ARGS[1]
    steady=solve_three_steady_states(default_parameters();verbose=false)
    expansion=stage3_expansion_context(steady)
    boom_x=reference_vector(dir,"original_boom_candidate")
    length(boom_x)==2*length(expansion.economy.B)-1 ||
        error("Original boom x has wrong dimension")
    # Warm-start expansion at the author's ACTUAL saved MATLAB boom solution.
    # Only a converged preliminary equilibrium may supply inherited states.
    boom=solve_stage3_flexible(expansion;x0=boom_x,
                                tolerance=1e-6,maxiter=5,verbose=false)
    println("PRE-SHOCK JULIA GE accepted residual=",
            boom.maximum_scaled_residual)
    source_mu=reference_matrix(dir,"shock_pre_distribution")
    metadata=reference_vector(dir,"shock_inherited_metadata")
    length(metadata)==7 || error("Expected seven source inherited scalars")
    failures=String[]
    for (name,lev,isfixed) in (("unleveraged",0.0,false),
                               ("leveraged",1.0/16.0,false),
                               ("fixed",1.0/16.0,true))
        ctx=stage3_reversal_context(steady,boom;
                                    shock_start=41,leverage_hat=lev)
        size(source_mu)==size(ctx.mu_pre) ||
            error("Author shock-date distribution differs in dimension")
        println("SOURCE INHERITANCE ",name,
                " max_initial_distribution_error=",
                maximum(abs.(ctx.mu_pre.-source_mu)),
                " Kdiff=",ctx.economy.initial.K-metadata[1],
                " qprev_diff=",ctx.economy.q_initial-metadata[2],
                " Bprev_diff=",ctx.economy.B_initial-metadata[3],
                " terms_prev_diff=",ctx.terms_pre-metadata[4])
        # Evaluate the author's exact pre-shock distribution, thus
        # isolating household-map differences from pre-shock distribution
        # numerical discretization. Do not overwrite production state.
        matched=Stage3Context(ctx.economy,copy(source_mu),
                              ctx.terms_pre,ctx.leverage_hat,
                              ctx.theta_nr,ctx.foreign,ctx.label)
        ok=check_author_solved_path(dir,name,matched;fixed=isfixed)
        ok || push!(failures,name)
    end
    isempty(failures) ||
        error("Julia does NOT reproduce stored author GE solution(s): " *
              join(failures,", ") *
              ". This is a numerical failure, not a nonexistent solution.")
    println("ALL THREE ACTUAL MATLAB SOLVED PATHS reproduced by Julia residual maps.")
end
main()
