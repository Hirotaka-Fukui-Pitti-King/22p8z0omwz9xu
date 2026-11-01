"""
Calibrate beta so that the household stationary net-asset demand equals
reference capital plus reference net foreign assets. Beta=0.9832 is the
author's warm-start benchmark, not an imposed result.
"""
function calibrate_beta(p::ModelParams{T},ref::ReferenceCalibration{T},
                        chain::ProductivityChain{T},gr::AssetGrids{T};
                        lo::T=T(0.960), hi::T=T(0.995)) where {T}
    B = ref.nfa_y*ref.Y
    aggregate = aggregate_prices(p,ref,ref.r,B)
    prices = HouseholdPrices(ref.r,aggregate.wage,aggregate.terms,
                             aggregate.L,aggregate.profits)
    function residual(beta::T)
        hh = solve_household(p,gr,chain,prices,beta)
        mu = stationary_distribution(p,gr,chain,hh).distribution
        assets = sum(mu .* hh.savings_fine)
        return (assets-(ref.K+B))/max(one(T),abs(ref.K+B))
    end
    result = bracketed_root(residual,lo,hi;tol=p.root_tol,maxiter=p.root_maxiter)
    return (beta=result.root,scaled_residual=result.residual,iterations=result.iterations)
end

"""
Solve the full stationary GE in r, nested around EGM and distribution.
Every root evaluation recomputes prices, household policies, and household
distribution for the given B. No stationary objects are hard-coded.
"""
function solve_steady_state(p::ModelParams{T},ref::ReferenceCalibration{T},
                            chain::ProductivityChain{T},gr::AssetGrids{T},
                            beta::T,B::T;
                            name::String="steady",lo::T=T(-0.005),
                            hi::T=T(0.055)) where {T}
    function residual(r::T)
        st = stationary_candidate(p,ref,chain,gr,beta,r,B;name=name)
        return st.asset_market_error/max(one(T),abs(st.K+B))
    end
    solution = bracketed_root(residual,lo,hi;tol=p.root_tol,maxiter=p.root_maxiter)
    st = stationary_candidate(p,ref,chain,gr,beta,solution.root,B;name=name)
    st_err = abs(st.asset_market_error)/max(one(T),abs(st.K+B))
    st_err <= p.root_tol ||
        throw(ErrorException("post-solve GE asset residual too high: $(st_err)"))
    return (state=st,root_iterations=solution.iterations)
end

function solve_three_steady_states(p::ModelParams{T}=default_parameters()) where {T}
    ref = baseline_reference(p)
    chain = rouwenhorst(p)
    grids = make_asset_grids(p,ref.K)
    calibration = calibrate_beta(p,ref,chain,grids)
    beta = calibration.beta
    B_cal = ref.nfa_y*ref.Y
    factor = (one(T)-p.credit_persistence^p.transition_horizon)/
             (one(T)-p.credit_persistence)
    B_final = -p.current_account_initial*ref.Y*factor
    calibrated = solve_steady_state(p,ref,chain,grids,beta,B_cal;name="calibrated")
    initial = solve_steady_state(p,ref,chain,grids,beta,zero(T);name="initial")
    final = solve_steady_state(p,ref,chain,grids,beta,B_final;name="final")
    return (parameters=p,reference=ref,chain=chain,grids=grids,
            beta_calibration=calibration,
            calibrated=calibrated.state,initial=initial.state,final=final.state,
            root_iterations=(calibrated=calibrated.root_iterations,
                             initial=initial.root_iterations,final=final.root_iterations))
end
