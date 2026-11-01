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
    result = solve_flexible_transition(inputs;
                                      tolerance=1e-6,maxiter=70,
                                      refresh_every=15,verbose=true)
    diagnostic = transition_diagnostics(result;compute_euler=true)
    outdir=normpath(joinpath(@__DIR__,"..","results"))
    mkpath(outdir)
    csv = export_transition_csv(result,joinpath(outdir,"stage2_transition_expansion.csv"))
    figs = write_transition_svg(result,outdir)
    diagnostic_file=joinpath(outdir,"stage2_equilibrium_diagnostics.csv")
    open(diagnostic_file,"w") do io
        println(io,"metric,value")
        for key in keys(diagnostic)
            println(io,string(key),",",getproperty(diagnostic,key))
        end
    end
    logfile=joinpath(outdir,"stage2_solver_log.txt")
    open(logfile,"w") do io
        println(io,"Stage 2 nonlinear solver")
        println(io,"Julia version: ",VERSION)
        println(io,"Project: ",Base.active_project())
        println(io,"Converged: ",result.converged)
        println(io,"Iterations: ",result.iterations)
        println(io,"Residual evaluations: ",result.evaluations)
        println(io,"Maximum scaled solver residual: ",result.maximum_scaled_residual)
        println(io,"Author labor exponent discrepancy remains: true")
        for key in keys(diagnostic)
            println(io,string(key)," = ",getproperty(diagnostic,key))
        end
    end
    for (key,value) in pairs(diagnostic)
        println("DIAGNOSTIC ",key," = ",value)
    end
    println("Saved ",csv)
    println("Saved ",diagnostic_file)
    println("Saved ",logfile)
    println("SVG: ",join(figs,", "))
    # This script never declares a numerical/paper replication pass on
    # solver convergence alone; values must meet independent conditions.
    if diagnostic.max_asset_market_level > 1e-4 ||
       diagnostic.max_labor_market_level > 1e-4 ||
       diagnostic.max_goods_market_level > 1e-4 ||
       diagnostic.max_probability_mass_deviation > 1e-8
        error("Independent GE residual test failed; inspect the exported diagnostics")
    end
    println("STAGE2: numerical equilibrium residuals passed preliminary level checks; compare the actual SVG/CSV with author Figures 1-3")
end

main()
