using FMR2020

function main()
    println("FMR-2020 Stage 2: nonlinear 200-quarter foreign credit expansion")
    println("Julia version: ",VERSION)
    println("Project: ",Base.active_project())
    println("Stage 1 calibration and three stationary equilibria:")
    steady = solve_three_steady_states(default_parameters();verbose=true)
    inputs = transition_inputs(steady)
    println("Transition periods: ",length(inputs.B),
            " | dimension: ",2*length(inputs.B)-1,
            " | initial B: ",inputs.B_initial,
            " | terminal B: ",inputs.B[end])
    println("Author original algorithm: aggregate K and L root, backward EGM, forward distribution")
    println("Author labor-supply exponent discrepancy: see docs/STAGE2_AUDIT.md")
    outdir=normpath(joinpath(@__DIR__,"..","results"))
    mkpath(outdir)
    result = try
        solve_flexible_transition(inputs;
                                  tolerance=1e-6,maxiter=70,
                                  refresh_every=15,verbose=true)
    catch err
        path=joinpath(outdir,"stage2_solver_log.txt")
        open(path,"w") do io
            println(io,"Stage 2 nonlinear equilibrium FAILED")
            println(io,"Julia: ",VERSION)
            println(io,"Error: ",sprint(showerror,err))
            println(io,"Do not report an unconverged path or plot guessed IRFs.")
        end
        println("Wrote failure log: ",path)
        rethrow()
    end
    diagnostic = transition_diagnostics(result;compute_euler=true)
    diagnostic_file=joinpath(outdir,"stage2_equilibrium_diagnostics.csv")
    open(diagnostic_file,"w") do io
        println(io,"metric,value")
        for key in keys(diagnostic)
            println(io,string(key),",",getproperty(diagnostic,key))
        end
    end
    for (key,value) in pairs(diagnostic)
        println("DIAGNOSTIC ",key," = ",value)
    end
    # A root return flag and the 399 residuals do not substitute for
    # independent household-aggregated INTERIOR goods clearing.
    # The original authors do not impose a terminal-period (T-th) goods
    # market condition: report that value separately, never hide it.
    valid = diagnostic.max_asset_market_level <= 1e-4 &&
            diagnostic.max_labor_market_level <= 1e-4 &&
            diagnostic.max_goods_market_level <= 1e-4 &&
            diagnostic.max_probability_mass_deviation <= 1e-8 &&
            diagnostic.terminal_capital_gap/
                max(1.0,abs(inputs.terminal.K)) <= 0.02 &&
            diagnostic.terminal_labor_gap/
                max(1.0,abs(inputs.terminal.L)) <= 0.02
    logpath=joinpath(outdir,"stage2_solver_log.txt")
    open(logpath,"w") do io
        println(io,"Stage 2 nonlinear solver")
        println(io,"Julia version: ",VERSION)
        println(io,"Project: ",Base.active_project())
        println(io,"Root converged: ",result.converged)
        println(io,"Independent interior GE accepted: ",valid)
        println(io,"Iterations: ",result.iterations)
        println(io,"Residual evaluations: ",result.evaluations)
        println(io,"Maximum scaled 399-equation residual: ",result.maximum_scaled_residual)
        println(io,"Authors' terminal consumption and wealth are boundary approximations, not extra T-th market equations.")
        println(io,"The terminal goods/budget gaps are REPORTED but not used in the 1:T-1 goods market gate.")
        for key in keys(diagnostic)
            println(io,string(key)," = ",getproperty(diagnostic,key))
        end
    end
    println("Saved ",diagnostic_file)
    println("Saved ",logpath)
    if !valid
        error("Independent INTERIOR general equilibrium residual test failed; " *
              "no new transition CSV or SVG was generated. See stage2 diagnostics.")
    end
    csv=export_transition_csv(result,joinpath(outdir,"stage2_transition_expansion.csv"))
    figs=write_transition_svg(result,outdir)
    println("Saved ",csv)
    println("SVG: ",join(figs,", "))
    println("STAGE2: numerical equilibrium residuals passed preliminary level checks; compare the actual SVG/CSV with author Figures 1-3")
end

main()
