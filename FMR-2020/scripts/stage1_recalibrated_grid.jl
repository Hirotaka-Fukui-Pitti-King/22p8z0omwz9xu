using FMR2020

"""
Recompute beta and ALL three r equilibria on each grid. This is not the
fixed-beta/price test in scripts/grid_refinement.jl; every grid independently
clears the asset market. These runs can be considerably more expensive.
"""
function run_grid(n_coarse::Int,n_fine::Int)
    println("START recalibrated GE grid ",n_coarse,"/",n_fine)
    p = default_parameters(n_coarse=n_coarse,n_fine=n_fine)
    result = solve_three_steady_states(p;verbose=true)
    rows = Tuple[]
    for st in (result.initial,result.calibrated,result.final)
        diag = steady_state_diagnostics(p,result.chain,result.grids,st)
        println("GRID ",n_coarse,"/",n_fine," ",st.name,
                " beta=",st.beta," r=",st.r,
                " K=",st.K," L=",st.L,
                " max Euler=",diag.household.max_euler_unconstrained,
                " asset scaled residual=",diag.scaled_asset_error)
        push!(rows,(n_coarse,n_fine,st.name,st.beta,st.r,st.B,st.K,st.L,
                    st.Y,st.consumption_home,st.aggregate_assets,
                    diag.scaled_asset_error,diag.scaled_goods_error,
                    diag.household.max_euler_unconstrained,
                    diag.distribution.invariance_error))
    end
    return rows
end

function main()
    dest = normpath(joinpath(@__DIR__,"..","results"))
    mkpath(dest)
    pairs = get(ENV,"FMR_LARGE_GRID","0") == "1" ?
        ((97,385),(193,769),(385,1537)) : ((97,385),(193,769))
    path = joinpath(dest,"stage1_recalibrated_grid.csv")
    open(path,"w") do io
        println(io,"n_coarse,n_fine,steady_state,beta,r,B,K,L,Y,C_H,A,scaled_asset_error,scaled_goods_error,max_coarse_euler,distribution_invariance")
        for (nc,nf) in pairs
            for row in run_grid(nc,nf)
                println(io,join(row,","))
            end
            flush(io)
        end
    end
    println("Saved ",path)
    println("Compare across grids within each steady_state; beta is recalibrated per grid.")
end

main()
