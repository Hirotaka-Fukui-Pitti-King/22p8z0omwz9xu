# Compare the ORIGINAL author's solved fixed-rate 599 vector to Julia
# prices, distribution, realized savings and policy arrays at identical
# K, L, terms, inherited wealth and shock distribution. NO root-finding.
using FMR2020
using LinearAlgebra
include(joinpath(@__DIR__,"evaluate_author_explicit_goods.jl"))

function finite_source_matrix(dir::String,label::String)
    file=joinpath(dir,label*".csv")
    isfile(file) || error("Missing exact MATLAB source audit fixture: "*file)
    lines=readlines(file)
    isempty(lines) && error("Empty MATLAB source audit fixture: "*file)
    nc=length(split(strip(first(lines)),","))
    result=Matrix{Float64}(undef,length(lines),nc)
    for (i,line) in enumerate(lines)
        cells=split(strip(line),",")
        length(cells)==nc || error("Ragged source audit fixture: "*file)
        for j in eachindex(cells)
            result[i,j]=parse(Float64,strip(cells[j]))
        end
    end
    return result
end

function compare_block(name::String,source::AbstractArray,
                       production::AbstractArray)
    size(source)==size(production) ||
        error("$name MATLAB $(size(source)) vs Julia $(size(production))")
    smask=isfinite.(source)
    jmask=isfinite.(production)
    nbad=count(smask .⊻ jmask)
    valid=smask .& jmask
    any(valid) || error("$name has no finite source-and-Julia comparison")
    gap=abs.(source[valid].-production[valid])
    maxgap=maximum(gap)
    idx=argmax(gap)
    println("FIXED MATCH ",name," finite_n=",count(valid),
            " nonfinite_mask_mismatch=",nbad,
            " max_abs_error=",maxgap,
            " finite_argmax=",idx)
    return maxgap
end

function main_fixed()
    length(ARGS)==1 ||
        error("Usage: julia --project=. scripts/diagnose_author_fixed_map.jl FIXTURE_DIR")
    dir=ARGS[1]
    st=solve_three_steady_states(default_parameters();verbose=false)
    mu=source_matrix(dir,"shock_pre_distribution")
    meta=source_vector(dir,"shock_inherited_metadata")
    ctx=source_reversal_context(st,mu,meta,1.0/16.0)
    x=source_vector(dir,"fixed_candidate")
    author_market=source_vector(dir,"fixed_source_residual")
    length(x)==599 && length(author_market)==599 ||
        throw(DimensionMismatch("Original fixed 599 author vectors required"))
    w=FMR2020.TransitionWorkspace(ctx.economy)
    res=zeros(Float64,599)
    FMR2020.stage3_fixed_residual!(res,x,ctx,w)
    raw=similar(res)
    stage3_unscale_residual!(raw,res,ctx;fixed=true)
    delta=abs.(raw.-author_market)
    idx=argmax(delta)
    println("ORIGINAL FIXED 599 SOURCE discrepancy max_raw=",delta[idx],
            " residual_index=",idx," source=",author_market[idx],
            " julia=",raw[idx])
    compare_block("interest rate",source_vector(dir,"fixed_src_rate"),w.rate)
    compare_block("wage",source_vector(dir,"fixed_src_wage"),w.wage)
    compare_block("firm profits",source_vector(dir,"fixed_src_profits"),w.profits)
    compare_block("capital q",source_vector(dir,"fixed_src_q"),w.q)
    compare_block("shock-date distribution",
                  finite_source_matrix(dir,"fixed_src_initial_distribution"),
                  view(w.distribution,:,:,1))
    compare_block("household aggregate savings",
                  source_vector(dir,"fixed_src_household_assets"),
                  view(w.aggregate_savings,1:199))
    for t in (1,9)
        compare_block("fine consumption at t=$t",
                      finite_source_matrix(dir,"fixed_src_consumption_t$t"),
                      view(w.fine_consumption,:,:,t))
        compare_block("fine savings at t=$t",
                      finite_source_matrix(dir,"fixed_src_savings_t$t"),
                      view(w.fine_savings,:,:,t))
    end
    println("FIXED SOURCE BLOCK AUDIT COMPLETE: all results are" *
            " SAME-PATH comparisons, NOT a claim of Julia equilibrium.")
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main_fixed()
end
