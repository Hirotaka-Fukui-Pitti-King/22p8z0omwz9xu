using FMR2020

function asset_supply_for_grid(ncoarse::Int,nfine::Int)
    p=default_parameters(n_coarse=ncoarse,n_fine=nfine,
                         egm_tol=1e-10,dist_tol=1e-10)
    reference=baseline_reference(p)
    process=rouwenhorst(p)
    grid=make_asset_grids(p,reference.K)
    B=reference.nfa_y*reference.Y
    aggregate=FMR2020.aggregate_prices(p,reference,reference.r,B)
    prices=HouseholdPrices(reference.r,aggregate.wage,aggregate.terms,
                           aggregate.L,aggregate.profits)
    hh=solve_household(p,grid,process,prices,0.9832)
    dist=stationary_distribution(p,grid,process,hh)
    assets=sum(dist.distribution .* hh.savings_fine)
    e=household_diagnostics(p,grid,process,prices,0.9832,hh)
    inv=distribution_invariance(dist.distribution,hh,grid,process)
    return (ncoarse,nfine,assets,sum(dist.distribution .* hh.consumption_fine),
            e.max_euler_unconstrained,e.max_kkt_violation,
            inv.invariance_error,hh.iterations,dist.iterations)
end

function main()
    sizes=((61,241),(97,385),(193,769))
    dest=joinpath(@__DIR__,"..","results")
    mkpath(dest)
    rows=map(pair->asset_supply_for_grid(pair[1],pair[2]),sizes)
    output=joinpath(dest,"stage1_grid_refinement.csv")
    open(output,"w") do io
        println(io,"n_coarse,n_fine,asset_supply,consumption_home,euler_error,kkt_error,invariance_error,egm_iterations,distribution_iterations")
        for row in rows
            println(io,join(row,","))
        end
    end
    for row in rows
        println("GRID ",row[1],"/",row[2]," assets=",row[3],
                " consumption=",row[4]," Euler=",row[5],
                " invariant=",row[7])
    end
    println("Difference large - middle (assets): ",rows[end][3]-rows[end-1][3])
    println("Difference middle - small (assets): ",rows[end-1][3]-rows[1][3])
    println("Saved ",output)
    println("These are FIXED-beta/price refinements, not re-calibrated GE comparisons.")
end

main()
