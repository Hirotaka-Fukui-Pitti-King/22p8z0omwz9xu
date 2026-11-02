# Compare one-step ORIGINAL-MATLAB Rouwenhorst, asset grids and EGM
# against source_matlab_egm_step on identical prescribed synthetic inputs.
# This is an independent MATLAB-Julia equality test, NOT a GE solve.
# Usage: julia --project=. scripts/compare_matlab_egm_fixture.jl path/to/fixture
using FMR2020

function read_numeric_csv(path::AbstractString)
    isfile(path) || error("MATLAB fixture missing: $path")
    rows=readlines(path)
    isempty(rows) && error("Empty MATLAB fixture: $path")
    cells=[split(strip(line),",") for line in rows]
    n=length(cells[1])
    n>0 && all(length(x)==n for x in cells) ||
        throw(DimensionMismatch("Inconsistent MATLAB fixture columns: $path"))
    mat=Matrix{Float64}(undef,length(rows),n)
    for (i,row) in enumerate(cells),j in 1:n
        cell=strip(row[j])
        mat[i,j]=uppercase(cell)=="NA" ? NaN : parse(Float64,cell)
    end
    return mat
end

function assert_fixture(label::String, ours::AbstractArray, theirs::AbstractArray;
                        atol::Float64=1e-8)
    size(ours)==size(theirs) ||
        error("$label dimensions differ: Julia=$(size(ours)) MATLAB=$(size(theirs))")
    omiss=isnan.(ours)
    mmiss=isnan.(theirs)
    omission_diff=count(omiss .!= mmiss)
    println(label,": undefined MATLAB=",count(mmiss),
            " undefined Julia=",count(omiss),
            " unmatched_undefined=",omission_diff)
    if omission_diff!=0
        println("MATLAB undefined indices: ",findall(mmiss))
        println("Julia undefined indices: ",findall(omiss))
        println("NaN disagreement indices: ",findall(omiss .!= mmiss))
        error("SOURCE PCHIP DOMAIN DISAGREEMENT in $label; " *
              "NaN status must match at every asset/productivity index")
    end
    valid=.!mmiss
    all(isfinite.(ours[valid])) && all(isfinite.(theirs[valid])) ||
        error("Nonfinite source-policy values outside expected NaN mask")
    any(valid) || error("No numerically valid source comparison nodes: $label")
    dif=zeros(Float64,size(ours))
    dif[valid] .=abs.(ours[valid].-theirs[valid])
    k=argmax(dif)
    gap=dif[k]
    println(label,": max_abs_error=",gap," index=",Tuple(k),
            " Julia=",ours[k]," MATLAB=",theirs[k])
    gap<=atol ||
        error("ORIGINAL-SOURCE NUMERIC DISAGREEMENT in $label: $gap > $atol")
end

function main()
    length(ARGS)==1 ||
        error("Pass MATLAB fixture folder from export_source_egm_fixture.m")
    dir=ARGS[1]
    params=default_parameters()
    ref=baseline_reference(params)
    gr=make_asset_grids(params,ref.K)
    ch=rouwenhorst(params)
    assert_fixture("original coarse asset nodes",reshape(gr.coarse_full,:,1),
                   read_numeric_csv(joinpath(dir,"coarse_grid.csv")))
    assert_fixture("original fine asset nodes",reshape(gr.fine_full,:,1),
                   read_numeric_csv(joinpath(dir,"fine_grid.csv")))
    assert_fixture("original productivity",reshape(ch.levels,:,1),
                   read_numeric_csv(joinpath(dir,"productivity.csv")))
    assert_fixture("original Markov probabilities",ch.transition,
                   read_numeric_csv(joinpath(dir,"transition.csv"));atol=1e-12)
    # Same positive increasing source future policy as the MATLAB fixture.
    # Keeps the tested EGM interpolation inside the endogenous-knot range.
    future=repeat(reshape(1.0 .+ 0.02 .* max.(gr.coarse_full,0.0),:,1),
                  1,params.n_s)
    step=source_matlab_egm_step(gr,ch,params,future,0.98,0.01,0.01,
                                1.0,1.8,1.0,0.25)
    assert_fixture("source coarse consumption",step.coarse_consumption,
                   read_numeric_csv(joinpath(dir,"coarse_consumption.csv")))
    assert_fixture("source fine consumption",step.fine_consumption,
                   read_numeric_csv(joinpath(dir,"fine_consumption.csv")))
    assert_fixture("source fine savings",step.fine_savings,
                   read_numeric_csv(joinpath(dir,"fine_savings.csv")))
    assert_fixture("source fine budget",step.budget_gap,
                   read_numeric_csv(joinpath(dir,"fine_budget.csv")))

    # This second experiment checks all fixed-point updates, rather
    # than merely one backward EGM evaluation.
    st=source_matlab_stationary_egm(gr,ch,params,0.98,0.01,
                                    1.0,1.8,1.0,0.25;
                                    tolerance=1e-9,maxiter=5000)
    println("SOURCE STATIONARY iterations=",st.iterations,
            " max_update=",st.update_error)
    assert_fixture("source stationary coarse policy",
                   st.coarse_consumption,
                   read_numeric_csv(joinpath(dir,"stationary_coarse_consumption.csv"));
                   atol=1e-6)
    assert_fixture("source stationary fine policy",
                   st.fine_consumption,
                   read_numeric_csv(joinpath(dir,"stationary_fine_consumption.csv"));
                   atol=1e-6)
    println("SOURCE EGM NUMERICAL COMPARISON PASS. This alone does NOT verify GE.")
end

main()
