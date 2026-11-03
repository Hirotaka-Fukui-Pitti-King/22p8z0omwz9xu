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
    # For the exact source-side PCHIP boundary and NA mask, use
    # the MATLAB-exported original grid endpoints, rather than
    # Julia-generated coordinates that differ by <= 1e-11.
    matlab_coarse=vec(read_numeric_csv(joinpath(dir,"coarse_grid.csv")))
    matlab_fine=vec(read_numeric_csv(joinpath(dir,"fine_grid.csv")))
    source_grid=AssetGrids(matlab_coarse,matlab_fine,gr.coarse,
                           matlab_fine,matlab_coarse[33],matlab_coarse[end])
    # Same positive increasing source future policy as the MATLAB fixture.
    # Keeps the tested EGM interpolation inside the endogenous-knot range.
    future=repeat(reshape(1.0 .+ 0.02 .* max.(gr.coarse_full,0.0),:,1),
                  1,params.n_s)
    step=source_matlab_egm_step(source_grid,ch,params,future,0.98,0.01,0.01,
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
    st=source_matlab_stationary_egm(source_grid,ch,params,0.98,0.01,
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
    # Original-source evaluation at the ACTUAL calibrated terminal price
    # vector. The synthetic test above alone cannot show how much the
    # existing 161-node Julia production policy differs at author prices.
    calibrated_input=joinpath(dir,"calibrated_input")
    calibrated_output=joinpath(dir,"calibrated_output")
    if isdir(calibrated_input) && isdir(calibrated_output)
        pv=vec(read_numeric_csv(joinpath(calibrated_input,"prices.csv")))
        length(pv)==9 || error("calibrated source-price vector must have 9 scalars")
        full_coarse=vec(read_numeric_csv(joinpath(calibrated_input,"coarse_grid.csv")))
        full_fine=vec(read_numeric_csv(joinpath(calibrated_input,"fine_grid.csv")))
        source_grid=AssetGrids(full_coarse,full_fine,gr.coarse,
                               full_fine,full_coarse[33],full_coarse[end])
        supplied_ex=vec(read_numeric_csv(joinpath(calibrated_input,"productivity.csv")))
        supplied_P=read_numeric_csv(joinpath(calibrated_input,"transition.csv"))
        supplied_chain=ProductivityChain(supplied_ex,supplied_P,ch.invariant)
        stc=source_matlab_stationary_egm(
            source_grid,supplied_chain,params,pv[1],pv[2],pv[3],
            pv[4],pv[5],pv[6];tolerance=1e-9,maxiter=5000)
        println("CALIBRATED SOURCE TERMINAL iterations=",stc.iterations,
                " update_error=",stc.update_error)
        assert_fixture("calibrated source terminal full coarse",
            stc.coarse_consumption,
            read_numeric_csv(joinpath(calibrated_output,"source_coarse_policy.csv"));
            atol=1e-6)
        assert_fixture("calibrated source terminal full fine",
            stc.fine_consumption,
            read_numeric_csv(joinpath(calibrated_output,"source_fine_policy.csv"));
            atol=1e-6)
        old_coarse=read_numeric_csv(joinpath(calibrated_input,
                                            "julia_active_coarse_policy.csv"))
        original_active=stc.coarse_consumption[33:end,:]
        size(original_active)==size(old_coarse) ||
            throw(DimensionMismatch("161-node terminal policy dimensions differ"))
        difference=abs.(original_active.-old_coarse)
        maxgap=maximum(difference)
        imax=argmax(difference)
        println("CALIBRATED TERMINAL 193-source versus 161-production " *
                "max_abs_difference=",maxgap,
                " active_index=",Tuple(imax),
                " source=",original_active[imax],
                " current=",old_coarse[imax])
        isfinite(maxgap) || error("Nonfinite source/production active comparison")
        old_fine=read_numeric_csv(joinpath(calibrated_input,
                                          "julia_active_fine_policy.csv"))
        original_fine=stc.fine_consumption[129:end,:]
        size(original_fine)==size(old_fine) ||
            throw(DimensionMismatch("641-node terminal fine policy dimensions differ"))
        vmask=isfinite.(original_fine)
        any(vmask) || error("No defined source terminal fine policy nodes")
        finitegap=maximum(abs.(original_fine[vmask].-old_fine[vmask]))
        println("CALIBRATED TERMINAL FINE source vs production finite_max_difference=",
                finitegap," source_nan_count=",count(.!vmask))
        isfinite(finitegap) || error("Nonfinite source/production fine comparison")
    else
        error("Calibrated terminal source fixture missing: expected both " *
              "calibrated_input and calibrated_output")
    end
    # Compare ORIGINAL source basefun.m revaluation of the REAL shock-date
    # 769x7 household distribution, for both roots and leverage regimes.
    fx_input=joinpath(dir,"fx_input")
    fx_output=joinpath(dir,"fx_output")
    isdir(fx_input) && isdir(fx_output) ||
        error("FX distribution source comparison fixture missing")
    fine=vec(read_numeric_csv(joinpath(fx_input,"fine_grid.csv")))
    mu_pre=read_numeric_csv(joinpath(fx_input,"mu_pre.csv"))
    metadata=read_numeric_csv(joinpath(fx_input,"metadata.csv"))
    matlab_fsolve=vec(read_numeric_csv(joinpath(fx_output,
                                               "original_fsolve_root.csv")))
    size(metadata)==(2,11) && size(mu_pre)==(length(fine),7) ||
        throw(DimensionMismatch("FX fixture metadata or distribution dimensions"))
    for i in 1:2
        dat=view(metadata,i,:)
        labels=i==1 ? "unleveraged" : "leveraged"
        f(x)=dat[6]+dat[7]/x-dat[8]/x^dat[9]
        abs(f(matlab_fsolve[i]))<1e-8 ||
            error("Original source fsolve trade root residual exceeded 1e-8")
        roots=FMR2020.foreign_debt_trade_roots(dat[6],dat[7],dat[8],dat[9])
        selected=abs(matlab_fsolve[i]-roots.low) <
                 abs(matlab_fsolve[i]-roots.high) ? "low" : "high"
        println("ORIGINAL MATLAB TRADE ROOT leverage=",dat[5],
                " from_x0_0.8=",matlab_fsolve[i],
                " selected_branch=",selected,
                " Julia_low=",roots.low," Julia_high=",roots.high)
        for (j,label) in enumerate(("low","high"))
            current=j==1 ? roots.low : roots.high
            result=zeros(Float64,size(mu_pre))
            FMR2020.revalue_foreign_debt!(result,mu_pre,fine,
                dat[1],dat[2],dat[3],dat[4],current,dat[5])
            author=read_numeric_csv(joinpath(fx_output,labels*"_"*
                                               label*"_distribution.csv"))
            assert_fixture("FX source basefun "*labels*" "*label,
                           result,author;atol=1e-9)
            abs(sum(result)-sum(mu_pre))<1e-10 ||
                error("FX source mass conservation failed")
        end
    end
    println("FX ORIGINAL BASEFUN VS JULIA NUMERIC COMPARISON PASS")
    # A separately labeled source-PCHIP-tail-extrapolated 399 MATLAB function is independently
    # executed once using the exported same candidate path, G0, B, and
    # calibrated 193-row terminal policy. It can be numerically
    # undefined; this MUST fail rather than being called a match.
    ge_output=joinpath(dir,"original_399_output")
    isdir(ge_output) ||
        error("Exact original MATLAB 399-residual output missing")
    assert_fixture("PCHIP-tail-extrapolated MATLAB raw 399 GE residual",
                   read_numeric_csv(joinpath(fx_input,"julia_raw_399.csv")),
                   read_numeric_csv(joinpath(ge_output,"matlab_raw_399.csv"));
                   atol=1e-6)
    for label in ("terms","rate","wage","aggregate_C")
        assert_fixture("tail-extrapolated source 399 "*label,
                       read_numeric_csv(joinpath(fx_input,"julia_"*label*".csv")),
                       read_numeric_csv(joinpath(ge_output,"matlab_"*label*".csv"));
                       atol=1e-6)
    end
    source_A=read_numeric_csv(joinpath(ge_output,"matlab_aggregate_A.csv"))
    # MATLAB's returned a_prime is target capital+bonds, NOT household
    # realized savings; fixture exporter reconstructs a_prime_t from
    # source asset-market residual R_K and the candidate K/B/q.
    assert_fixture("tail-extrapolated source 399 aggregate savings",
                   read_numeric_csv(joinpath(fx_input,"julia_aggregate_A.csv")),
                   source_A;atol=1e-6)
    println("DIRECT TAIL-EXTRAPOLATED MATLAB 399-EQUATION RESIDUAL COMPARISON PASS")
    println("SOURCE EGM NUMERICAL COMPARISON PASS. This alone does NOT verify GE.")
end

main()
