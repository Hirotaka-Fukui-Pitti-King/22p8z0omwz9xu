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
    nbad=count(xor.(smask,jmask))
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

function compare_terminal(name::String,source::Matrix{Float64},
                          julia_policy::AbstractMatrix{Float64})
    size(source,2)==size(julia_policy,2) ||
        throw(DimensionMismatch("$name productivity dimension"))
    if size(source,1)>size(julia_policy,1)
        # Original 193/769 grid includes 32/128 pre-borrowing knots;
        # Julia stationary household stores only 161/641 active knots.
        active=view(source,(end-size(julia_policy,1)+1):size(source,1),:)
        compare_block(name*" (active borrowing support)",active,julia_policy)
    elseif size(source,1)==size(julia_policy,1)
        compare_block(name,source,julia_policy)
    else
        throw(DimensionMismatch("$name original smaller than Julia policy"))
    end
end

function compare_occupied_policy(name::String,source::Matrix{Float64},
                                 production::AbstractMatrix{Float64},
                                 mu::AbstractMatrix{Float64})
    size(source)==size(production)==size(mu) ||
        throw(DimensionMismatch("$name occupied policy dimensions"))
    gap=zero(Float64)
    weighted=zero(Float64)
    mismatch_mass=zero(Float64)
    maxat=(0,0)
    @inbounds for j in axes(mu,2),i in axes(mu,1)
        mass=mu[i,j]
        mass==0 && continue
        if !(isfinite(source[i,j]) && isfinite(production[i,j]))
            mismatch_mass+=mass
            continue
        end
        difference=abs(source[i,j]-production[i,j])
        weighted+=mass*difference
        if difference>gap
            gap=difference
            maxat=(i,j)
        end
    end
    println("OCCUPIED FIXED ",name," max_abs_error=",gap,
            " max_at_grid_state=",maxat,
            " mass_weighted_abs_error=",weighted,
            " probability_with_undefined_policy=",mismatch_mass)
    return gap
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
    compare_terminal("terminal coarse cpol",source_matrix(dir,"fixed_src_terminal_coarse"),
                     view(w.coarse_consumption,:,:,size(w.coarse_consumption,3)))
    compare_terminal("terminal fine c_fine",source_matrix(dir,"fixed_src_terminal_fine"),
                     st.calibrated.household.consumption_fine)
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
    for t in (199,9,1)
        src_c=finite_source_matrix(dir,"fixed_src_consumption_t$t")
        src_a=finite_source_matrix(dir,"fixed_src_savings_t$t")
        actual_c=view(w.fine_consumption,:,:,t)
        actual_a=view(w.fine_savings,:,:,t)
        compare_block("fine consumption at t=$t",src_c,actual_c)
        compare_block("fine savings at t=$t",src_a,actual_a)
        compare_occupied_policy("consumption t=$t",src_c,actual_c,
                                view(w.distribution,:,:,t))
        compare_occupied_policy("savings t=$t",src_a,actual_a,
                                view(w.distribution,:,:,t))
    end
    println("FIXED SOURCE BLOCK AUDIT COMPLETE: all results are" *
            " SAME-PATH comparisons, NOT a claim of Julia equilibrium.")
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main_fixed()
end
