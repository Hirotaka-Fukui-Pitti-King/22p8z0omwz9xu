"""
    exact_static_levels(calibration, K, X, mu_previous)

Exact ex-post Nirei (2015) equations (7),(8),(11),(12) for a GIVEN
predetermined K and positive finite-N firm investment X.
Does NOT solve conditional pricing equations (9),(17), nor stochastic RE.
"""
function exact_static_levels(calibration::CalibrationResult, K::Real, X::Real,
                             mu_previous::Real; iterations::Int=110)
    p=calibration.parameters
    c=calibration.constants
    capital=Float64(K)
    invest=Float64(X)
    mu_prev=Float64(mu_previous)
    isfinite(capital) && capital>0.0 || throw(DomainError(K,"positive K required"))
    isfinite(invest) && invest>0.0 || throw(DomainError(X,"ES1 requires positive X"))
    isfinite(mu_prev) && mu_prev>0.0 || throw(DomainError(mu_previous,"positive previous mu required"))
    0.0<p.alpha<1.0 && p.zeta>1.0-p.alpha && p.psi>0.0 && c.c_L>0.0 ||
        throw(ArgumentError("monotonicity requires 0<alpha<1,zeta>1-alpha,psi,c_L>0"))
    iterations>=50 || throw(ArgumentError("iterations must be >=50"))
    Lmax=p.psi^(-inv(p.zeta))
    Ymax=capital^p.alpha*Lmax^(1.0-p.alpha)
    invest<Ymax || throw(DomainError(invest,"no admissible C>0 and leisure>0"))
    Lmin=(invest/capital^p.alpha)^(inv(1.0-p.alpha))
    0.0<Lmin<Lmax || throw(DomainError(Lmin,"invalid labor bracket"))
    # MRS and labor demand imply monotone investment as a function of L.
    implied_X(labor) = capital^p.alpha*labor^(1.0-p.alpha)*
        (1.0+c.c_L/p.zeta-c.c_L/(p.psi*p.zeta*labor^p.zeta))
    lo,hi=Lmin,Lmax
    for _ in 1:iterations
        mid=lo+(hi-lo)/2.0
        if implied_X(mid)<invest
            lo=mid
        else
            hi=mid
        end
    end
    labor=lo+(hi-lo)/2.0
    output=capital^p.alpha*labor^(1.0-p.alpha)
    consumption=output-invest
    leisure=1.0-p.psi*labor^p.zeta
    consumption>0.0 && leisure>0.0 ||
        throw(DomainError((consumption,leisure),"outside utility domain"))
    wage=c.c_L*output/labor
    mu=consumption^(-p.sigma)*leisure^(1.0-p.sigma)
    R=mu_prev/(p.beta*mu)
    all(v->isfinite(v)&&v>0.0,(output,consumption,labor,wage,mu,R)) ||
        throw(DomainError((output,consumption,labor,wage,mu,R),"nonpositive/nonfinite levels"))
    resources_residual=(output-consumption-invest)/output
    production_residual=output/(capital^p.alpha*labor^(1.0-p.alpha))-1.0
    labor_residual=wage*leisure/(p.psi*p.zeta*consumption*labor^(p.zeta-1.0))-1.0
    wage_residual=wage*labor/(c.c_L*output)-1.0
    pricing_kernel_residual=p.beta*mu*R/mu_prev-1.0
    maximum(abs,(resources_residual,production_residual,labor_residual,
                 wage_residual,pricing_kernel_residual))<1e-9 ||
        error("nonlinear static residual exceeds 1e-9")
    return (K=capital,X=invest,Y=output,C=consumption,L=labor,
            w=wage,mu=mu,R=R,resources_residual=resources_residual,
            production_residual=production_residual,
            labor_residual=labor_residual,wage_residual=wage_residual,
            pricing_kernel_residual=pricing_kernel_residual)
end

"""
    simulate_exact_level_period!(state)

Original 2015 3-choice finite-N capital and ES1; nonlinear exact level
household/market-clearing conditions and next-period realized mu.
Appendix A remains a local continuum forecast (NOT full stochastic RE).
"""
function simulate_exact_level_period!(state::SimulationState)
    ss=state.calibration.steady_state
    t=state.period
    prior_K=state.K
    prior_mu=state.mu_prior
    expected=compute_expected_equilibrium(state.re,state.calibration,prior_K,prior_mu)
    draw_productivity!(state.firms,state.rngs,state.calibration)
    eq=solve_equilibrium!(state.workspace,state.firms,state.calibration,
                           expected.X_next,prior_K)
    # Production at t depends on K_t, NOT selected K_(t+1).
    levels=exact_static_levels(state.calibration,prior_K,eq.X,prior_mu)
    stable_gap=log(levels.mu/ss.mu)-state.re.h*log(eq.K_next/ss.K)
    state.K=eq.K_next
    state.mu_prior=levels.mu
    state.period+=1
    return (t=t,K_prior=prior_K,K_next=eq.K_next,X=eq.X,
            expected_X=expected.X_next,epsilon=eq.epsilon,
            Y=levels.Y,C=levels.C,L=levels.L,w=levels.w,
            mu=levels.mu,R=levels.R,
            firm_capital_residual=eq.capital_residual,
            firm_investment_residual=eq.investment_residual,
            resources_residual=levels.resources_residual,
            labor_residual=levels.labor_residual,
            production_residual=levels.production_residual,
            pricing_kernel_residual=levels.pricing_kernel_residual,
            stable_eigenvector_gap=stable_gap)
end
