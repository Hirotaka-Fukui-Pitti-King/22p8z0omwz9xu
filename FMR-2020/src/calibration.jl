"""
Calibrate beta so that the household stationary net-asset demand equals
reference capital plus reference net foreign assets. Beta=0.9832 is the
author's warm-start benchmark, not an imposed result.
"""
function calibrate_beta(p::ModelParams{T},ref::ReferenceCalibration{T},
                        chain::ProductivityChain{T},gr::AssetGrids{T};
                        lo::T=T(0.980), hi::T=T(0.985)) where {T}
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
Find a sign-changing GE asset-market bracket *around the benchmark return*.
Unlike a fixed two-endpoint bisection, probing r close to the impatience bound
can push optimal savings beyond the finite distribution grid even when the
equilibrium at an interior rate is well-defined.

Only the explicitly identified upper-grid-domain failure is a reason to skip
an exploratory point. Other errors propagate. A valid sign-changing bracket
is required; this never fabricates a truncated stationary distribution.
"""
function equilibrium_rate_bracket(f, center::T, lo::T, hi::T;
                                  tolerance::T=T(1e-8),
                                  initial_step::T=T(0.0005),
                                  expansion::T=T(1.6),
                                  maximum_steps::Int=20) where {T<:AbstractFloat}
    lo <= center <= hi || throw(ArgumentError("reference rate must be in the search interval"))
    initial_step > zero(T) || throw(ArgumentError("initial_step must be positive"))
    expansion > one(T) || throw(ArgumentError("expansion must exceed 1"))
    middle = f(center)
    isfinite(middle) || throw(DomainError(middle,"nonfinite benchmark residual"))
    abs(middle) <= tolerance &&
        return (lower=center,upper=center,f_lower=middle,f_upper=middle,probes=1)
    samples = [(center,middle)]
    invalid_probes = 0
    for iteration in 1:maximum_steps
        shift = initial_step*expansion^(iteration-1)
        # Prefer the direction indicated by the usual increasing GE excess
        # asset-supply residual, but check both sides and do not assume
        # monotonicity without evaluating it.
        points = middle < 0 ? (min(hi,center+shift),max(lo,center-shift)) :
                              (max(lo,center-shift),min(hi,center+shift))
        for rate in points
            any(pair->pair[1] == rate,samples) && continue
            value = try
                f(rate)
            catch err
                if err isa DomainError &&
                   occursin("fine-grid policy exceeds distribution grid upper bound",
                            sprint(showerror,err))
                    invalid_probes += 1
                    continue
                end
                rethrow()
            end
            isfinite(value) || throw(DomainError(value,"nonfinite GE asset residual"))
            abs(value) <= tolerance &&
                return (lower=rate,upper=rate,f_lower=value,f_upper=value,
                        probes=length(samples)+1)
            push!(samples,(rate,value))
            sort!(samples,by=first)
            for k in 1:(length(samples)-1)
                left=samples[k]
                right=samples[k+1]
                if signbit(left[2]) != signbit(right[2])
                    return (lower=left[1],upper=right[1],
                            f_lower=left[2],f_upper=right[2],
                            probes=length(samples))
                end
            end
        end
        samples[1][1] <= lo && samples[end][1] >= hi && break
    end
    throw(ArgumentError(
        "No valid sign-changing GE interest-rate bracket: benchmark r=$(center), " *
        "F(center)=$(middle), explored r in [" *
        "$(samples[1][1]), $(samples[end][1])], skipped " *
        "$(invalid_probes) upper-grid-invalid candidates. " *
        "Check economics and asset upper tail before changing any bounds."))
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
    # A finite-asset stationary economy should not be evaluated at a return
    # above the natural impatience bound; otherwise the upper asset grid
    # can become an artificial destination for all savers.
    upper = min(hi,inv(beta)-one(T)-T(1e-4))
    lo < upper || throw(ArgumentError("no admissible positive interest-rate bracket"))
    function residual(r::T)
        st = stationary_candidate(p,ref,chain,gr,beta,r,B;name=name)
        return st.asset_market_error/max(one(T),abs(st.K+B))
    end
    bracket = equilibrium_rate_bracket(residual,ref.r,lo,upper;
                                       tolerance=p.root_tol)
    solution = if bracket.lower == bracket.upper
        (root=bracket.lower,residual=bracket.f_lower,iterations=0)
    else
        bracketed_root(residual,bracket.lower,bracket.upper;
                       tol=p.root_tol,maxiter=p.root_maxiter)
    end
    st = stationary_candidate(p,ref,chain,gr,beta,solution.root,B;name=name)
    st_err = abs(st.asset_market_error)/max(one(T),abs(st.K+B))
    st_err <= p.root_tol ||
        throw(ErrorException("post-solve GE asset residual too high: $(st_err)"))
    return (state=st,root_iterations=solution.iterations)
end

function solve_three_steady_states(p::ModelParams{T}=default_parameters();
                                   verbose::Bool=false) where {T}
    ref = baseline_reference(p)
    chain = rouwenhorst(p)
    grids = make_asset_grids(p,ref.K)
    verbose && println("STAGE1: beta calibration at benchmark r=",ref.r)
    calibration = calibrate_beta(p,ref,chain,grids)
    verbose && println("STAGE1: calibrated beta=",calibration.beta,
                       "; scaled asset residual=",calibration.scaled_residual)
    beta = calibration.beta
    B_cal = ref.nfa_y*ref.Y
    factor = (one(T)-p.credit_persistence^p.transition_horizon)/
             (one(T)-p.credit_persistence)
    B_final = -p.current_account_initial*ref.Y*factor
    verbose && println("STAGE1: calibrated steady-state GE (B=",B_cal,")")
    calibrated = solve_steady_state(p,ref,chain,grids,beta,B_cal;name="calibrated")
    verbose && println("STAGE1: calibrated r=",calibrated.state.r,
                       "; asset residual=",calibrated.state.asset_market_error)
    verbose && println("STAGE1: initial steady-state GE (B=0)")
    initial = solve_steady_state(p,ref,chain,grids,beta,zero(T);name="initial")
    verbose && println("STAGE1: initial r=",initial.state.r,
                       "; asset residual=",initial.state.asset_market_error)
    verbose && println("STAGE1: final steady-state GE (B=",B_final,")")
    final = solve_steady_state(p,ref,chain,grids,beta,B_final;name="final")
    verbose && println("STAGE1: final r=",final.state.r,
                       "; asset residual=",final.state.asset_market_error)
    return (parameters=p,reference=ref,chain=chain,grids=grids,
            beta_calibration=calibration,
            calibrated=calibrated.state,initial=initial.state,final=final.state,
            root_iterations=(calibrated=calibrated.root_iterations,
                             initial=initial.root_iterations,final=final.root_iterations))
end
