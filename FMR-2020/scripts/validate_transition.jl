using FMR2020

function main()
    println("FMR-2020 Stage 2: stationary-limit and live backward/forward smoke")
    println("Julia version: ",VERSION)
    p=default_parameters(n_coarse=61,n_fine=241,egm_tol=1e-9,dist_tol=1e-9)
    ref=baseline_reference(p)
    process=rouwenhorst(p)
    grids=make_asset_grids(p,ref.K)
    B=ref.nfa_y*ref.Y
    prices=FMR2020.aggregate_prices(p,ref,ref.r,B)
    household_prices=HouseholdPrices(ref.r,prices.wage,prices.terms,
                                     prices.L,prices.profits)
    hh=solve_household(p,grids,process,household_prices,0.9832)
    mu=stationary_distribution(p,grids,process,hh).distribution
    steady=FMR2020.stationary_candidate(p,ref,process,grids,0.9832,ref.r,B;
                                        name="benchmark-smoke")
    input=TransitionInputs(p,ref,process,grids,steady,steady,fill(B,8),
                           ones(Float64,8),17.0,1.0,B)
    check=transition_stationary_check(input)
    println("Stationary-limit residuals (including author labor formula mismatch): ",check)
    x=initial_transition_guess(input)
    @assert length(x)==15
    evaluated=evaluate_transition!(input,x)
    w=evaluated.workspace
    println("Shock-free transition max scaled residual: ",maximum(abs,evaluated.residual))
    println("Transition mass max error: ",maximum(w.distribution_error))
    println("Capital price max gap from 1: ",maximum(abs.(w.q .- 1)))
    println("Stationary benchmark labor: ",steady.L,
            "; implied labor: ",w.labor_implied[1])
    maximum(w.distribution_error)<1e-8 ||
        error("probability mass is not preserved")
    all(isfinite,evaluated.residual) || error("nonfinite transition residual")
    println("SMOKE PASS: prices, backward EGM, forward distribution and residual are executable.")
    println("IMPORTANT: SMOKE PASS is not a solved 200-quarter transition.")
end

main()
