using FMR2020

function main()
    println("FMR-2020 Stage 2: stationary-limit and live backward/forward smoke")
    println("Julia version: ",VERSION)
    p=default_parameters(n_coarse=61,n_fine=241,egm_tol=1e-9,dist_tol=1e-9)
    println("Solving TRUE stationary GE on the smoke-test grids")
    ge=solve_three_steady_states(p;verbose=true)
    # Calibrated state has labor normalization L=1. It is the clean
    # constant-path test for the author's labor closure; the uncalibrated
    # stationary_candidate(beta=0.9832) is NOT an asset-clearing GE.
    steady=ge.calibrated
    B=steady.B
    println("Stationary calibration asset error: ",
            steady.asset_market_error)
    input=TransitionInputs(p,ge.reference,ge.chain,ge.grids,steady,steady,
                           fill(B,8),ones(Float64,8),17.0,1.0,B)
    check=transition_stationary_check(input)
    println("Stationary-limit residuals (including author labor formula mismatch): ",check)
    x=initial_transition_guess(input)
    @assert length(x)==15
    evaluated=evaluate_transition!(input,x)
    w=evaluated.workspace
    # Test actual household-aggregated goods clearing on interior dates.
    # The source author does not include an independent T-th goods
    # condition; diagnose its boundary mismatch separately.
    constructed=FMR2020.TransitionSolveResult(input,w,copy(x),
                                              copy(evaluated.residual),0,
                                              w.function_evals,
                                              maximum(abs,evaluated.residual),true)
    goodscheck=transition_diagnostics(constructed;compute_euler=false)
    println("Shock-free max interior goods residual: ",
            goodscheck.max_goods_market_level,
            " at t=",goodscheck.max_goods_market_period)
    println("Shock-free terminal goods gap (NOT an interior GE equation): ",
            goodscheck.terminal_goods_gap_unconstrained)
    goodscheck.max_goods_market_level<1e-4 ||
        error("Constant-path interior goods-market clearing failed")
    println("Shock-free transition max scaled residual: ",maximum(abs,evaluated.residual))
    println("Transition mass max error: ",maximum(w.distribution_error))
    println("Capital price max gap from 1: ",maximum(abs.(w.q .- 1)))
    println("Stationary benchmark labor: ",steady.L,
            "; implied labor: ",w.labor_implied[1])
    maximum(w.distribution_error)<1e-8 ||
        error("probability mass is not preserved")
    maximum(abs,evaluated.residual)<1e-4 ||
        error("Solved stationary GE should satisfy the shock-free dynamic limit; " *
              "inspect capital/asset, labor, and terminal timing")
    all(isfinite,evaluated.residual) || error("nonfinite transition residual")
    println("SMOKE PASS: prices, backward EGM, forward distribution and residual are executable.")
    println("IMPORTANT: SMOKE PASS is not a solved 200-quarter transition.")
end

main()
