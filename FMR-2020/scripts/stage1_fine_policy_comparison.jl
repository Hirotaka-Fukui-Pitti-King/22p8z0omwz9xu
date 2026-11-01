using FMR2020

"""
Controlled comparison at the *same* calibrated benchmark beta and aggregate
prices. Only the fine-grid EGM policy evaluation changes. This isolates
off-grid policy accuracy; it does NOT recalibrate the three steady states.
"""
function compare_methods()
    p = default_parameters()
    reference = baseline_reference(p)
    chain = rouwenhorst(p)
    gr = make_asset_grids(p,reference.K)
    cal = calibrate_beta(p,reference,chain,gr)
    beta = cal.beta
    B = reference.nfa_y*reference.Y
    outdir = normpath(joinpath(@__DIR__,"..","results"))
    mkpath(outdir)
    output = joinpath(outdir,"stage1_fine_policy_comparison.csv")
    rows = Tuple[]
    for method in (:coarse_savings,:direct_endogenous)
        st = stationary_candidate(p,reference,chain,gr,beta,reference.r,B;
                                  name=string(method),fine_method=method)
        hd = steady_state_diagnostics(p,chain,gr,st)
        acc = euler_accuracy_map(p,gr,chain,st)
        s = acc.summary
        row = (string(method),beta,st.r,st.aggregate_assets,
               st.asset_market_error,hd.household.max_budget_error,
               hd.household.max_fine_budget_error,
               s.max_unconstrained,s.max_kkt_violation,
               s.probability_weighted_mean,s.probability_weighted_p99,
               s.mass_in_top_ten_percent_asset_grid,
               hd.distribution.invariance_error)
        push!(rows,row)
        println("METHOD ",method,": asset demand=",st.aggregate_assets,
                ", asset residual=",st.asset_market_error,
                ", max fine Euler=",s.max_unconstrained,
                ", weighted fine Euler=",s.probability_weighted_mean,
                ", mass p99=",s.probability_weighted_p99,
                ", fine budget=",hd.household.max_fine_budget_error)
    end
    open(output,"w") do io
        println(io,"fine_method,beta,r,aggregate_assets,asset_market_residual,max_coarse_budget,max_fine_budget,max_fine_euler,max_kkt,weighted_mean_euler,mass_p99_euler,upper_asset_tail_mass,distribution_invariance")
        for row in rows
            println(io,join(row,","))
        end
    end
    println("Saved ",output)
    println("Decision: compare Euler accuracy and asset-demand shifts BEFORE changing the default method.")
end

compare_methods()
