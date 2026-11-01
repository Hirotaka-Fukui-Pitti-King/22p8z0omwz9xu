using FMR2020

function main()
    p = default_parameters(n_coarse=61,n_fine=241,
                           egm_tol=1e-9,dist_tol=1e-9)
    reference = baseline_reference(p)
    chain = rouwenhorst(p)
    grids = make_asset_grids(p,reference.K)
    println("SMOKE productivity invariant mass: ",sum(chain.invariant))
    println("SMOKE normalized productivity: ",sum(chain.invariant .* chain.levels))
    println("SMOKE coarse/fine active sizes: ",
            length(grids.coarse)," / ",length(grids.fine))
    B = reference.nfa_y*reference.Y
    prices = FMR2020.aggregate_prices(p,reference,reference.r,B)
    household_prices = HouseholdPrices(reference.r,prices.wage,prices.terms,
                                       prices.L,prices.profits)
    household = solve_household(p,grids,chain,household_prices,0.9832)
    measure = stationary_distribution(p,grids,chain,household)
    hd = household_diagnostics(p,grids,chain,household_prices,0.9832,household)
    dd = distribution_invariance(measure.distribution,household,grids,chain)
    println("SMOKE EGM iterations: ",household.iterations)
    println("SMOKE EGM update residual: ",household.update_error)
    println("SMOKE Euler unconstrained residual: ",hd.max_euler_unconstrained)
    println("SMOKE Euler KKT residual: ",hd.max_kkt_violation)
    println("SMOKE budget error: ",hd.max_budget_error)
    println("SMOKE distribution iterations: ",measure.iterations)
    println("SMOKE distribution residual: ",dd.invariance_error)
    println("SMOKE mass: ",dd.mass)
    println("SMOKE aggregate assets: ",sum(measure.distribution .* household.savings_fine))
    abs(dd.mass-1) < 1e-8 || error("Smoke mass test failed")
    hd.max_budget_error < 1e-7 || error("Smoke budget test failed")
    println("SMOKE PASS (this is NOT a three-steady-state GE validation)")
end

main()
