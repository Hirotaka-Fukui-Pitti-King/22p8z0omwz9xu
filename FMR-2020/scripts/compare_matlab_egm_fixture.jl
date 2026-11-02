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
        mat[i,j]=parse(Float64,strip(row[j]))
    end
    return mat
end

function assert_fixture(label::String, ours::AbstractArray, theirs::AbstractArray;
                        atol::Float64=1e-8)
    size(ours)==size(theirs) ||
        error("$label dimensions differ: Julia=$(size(ours)) MATLAB=$(size(theirs))")
    dif=abs.(ours.-theirs)
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
    future=fill(10.0,193,params.n_s)
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
    println("SOURCE EGM NUMERICAL COMPARISON PASS. This alone does NOT verify GE.")
end

main()
