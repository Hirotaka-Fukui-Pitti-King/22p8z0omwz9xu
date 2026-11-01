using FMR2020

function summary_row(st::SteadyState)
    return (st.name,st.beta,st.r,st.B,st.K,st.L,st.Y,st.wage,
            st.terms,st.investment,st.capital_price,st.rental,st.profits,
            st.consumption_home,st.consumption_foreign,st.aggregate_assets,
            st.asset_market_error,st.goods_market_error,st.union_error,
            st.distribution_iterations,st.household.iterations)
end

function main()
    println("FMR-2020 Stage 1: three steady states")
    println("Julia version: ",VERSION)
    println("Project: ",Base.active_project())
    p = default_parameters()
    solution = solve_three_steady_states(p;verbose=true)
    println("Beta calibration: ",solution.beta_calibration)
    println("Reference nfa_y: ",solution.reference.nfa_y)
    println("Reference k_y: ",solution.reference.k_y)
    dest = joinpath(@__DIR__,"..","results")
    mkpath(dest)
    names = ("initial","calibrated","final")
    states = (solution.initial,solution.calibrated,solution.final)
    csv = joinpath(dest,"stage1_steady_states.csv")
    open(csv,"w") do io
        println(io,"name,beta,r,B,K,L,Y,w,x,I,q,rK,profits,C_H,C_F,A,asset_residual,goods_residual,union_residual,distribution_iterations,egm_iterations")
        for st in states
            println(io,join(summary_row(st),","))
        end
    end
    logpath = joinpath(dest,"stage1_diagnostics.txt")
    open(logpath,"w") do io
        println(io,"Julia version: ",VERSION)
        println(io,"Beta: ",solution.beta_calibration)
        println(io,"Root iterations: ",solution.root_iterations)
        for (name,st) in zip(names,states)
            diag = steady_state_diagnostics(p,solution.chain,solution.grids,st)
            println(io,"[$name] ",diag)
            println("[$name] asset=",diag.scaled_asset_error,
                    " goods=",diag.scaled_goods_error,
                    " Euler=",diag.household.max_euler_unconstrained,
                    " mass=",diag.distribution.mass)
            if diag.scaled_asset_error > p.root_tol
                error("Stationary asset market failed for $name")
            end
        end
    end
    println("Created ",csv)
    println("Created ",logpath)
    println("STAGE1 NUMERICAL RUN COMPLETE; inspect diagnostics before acceptance")
end

main()
