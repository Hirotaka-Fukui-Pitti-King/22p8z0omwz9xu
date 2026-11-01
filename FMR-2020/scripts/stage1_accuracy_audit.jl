using FMR2020

function main()
    println("FMR-2020 Stage 1 — state-by-state fine-grid Euler audit")
    println("Julia version: ",VERSION)
    solution = solve_three_steady_states(default_parameters();verbose=true)
    dest = normpath(joinpath(@__DIR__,"..","results"))
    mkpath(dest)
    path = joinpath(dest,"stage1_euler_by_state.csv")
    summary_path = joinpath(dest,"stage1_euler_summary.csv")
    open(path,"w") do io
        println(io,"steady_state,asset_index,productivity_index,asset,productivity,mass,next_asset,consumption,relative_euler_gap,euler_kkt_violation,binding")
        open(summary_path,"w") do so
            println(so,"steady_state,max_unconstrained,max_kkt_violation,probability_weighted_mean,probability_weighted_unconstrained,probability_weighted_p99,mass_in_top_ten_percent_asset_grid,max_state_probability,total_mass,n_active_states,n_grid_states")
            for st in (solution.initial,solution.calibrated,solution.final)
                diagnostic = euler_accuracy_map(solution.parameters,solution.grids,
                                                solution.chain,st)
                s = diagnostic.summary
                println(so,join((st.name,s.max_unconstrained,s.max_kkt_violation,
                                 s.probability_weighted_mean,
                                 s.probability_weighted_unconstrained,
                                 s.probability_weighted_p99,
                                 s.mass_in_top_ten_percent_asset_grid,
                                 s.max_state_probability,s.total_mass,
                                 s.n_active_states,s.n_grid_states),","))
                println("STATE ",st.name,": fine-grid max Euler=",s.max_unconstrained,
                        " weighted mean=",s.probability_weighted_mean,
                        " mass-p99=",s.probability_weighted_p99,
                        " mass at upper 10% asset grid=",
                        s.mass_in_top_ten_percent_asset_grid)
                for record in diagnostic.records
                    println(io,join((st.name,record.asset_index,
                                     record.productivity_index,record.asset,
                                     record.productivity,record.stationary_mass,
                                     record.next_asset,record.consumption,
                                     record.relative_euler_gap,
                                     record.euler_kkt_violation,record.binding),","))
                end
                for record in Iterators.take(diagnostic.ranked,5)
                    println("  largest residual: asset=",record.asset,
                            " productivity=",record.productivity,
                            " mass=",record.stationary_mass,
                            " next_asset=",record.next_asset,
                            " binding=",record.binding,
                            " violation=",record.euler_kkt_violation)
                end
            end
        end
    end
    println("Saved ",path)
    println("Saved ",summary_path)
    println("Diagnostic only: do not infer accuracy from market clearing alone.")
end

main()
