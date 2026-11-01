"""
Model parameters. All rates are quarterly and B is home-good-valued net foreign assets.
A negative B represents net foreign borrowing. Numerical defaults track the
published author code; they are not calibrated numerical results.
"""
Base.@kwdef struct ModelParams{T<:AbstractFloat}
    alpha::T = 0.33
    delta::T = 0.02
    chi::T = 0.60
    gamma::T = 1.0
    theta::T = 1.0
    theta_star::T = 3.0
    epsilon_f::T = 10.0
    epsilon_w::T = 10.0
    varphi::T = 0.50
    rho_s::T = 0.97
    sigma_log_stationary::T = 0.84
    n_s::Int = 7
    n_coarse::Int = 193
    n_fine::Int = 769
    grid_min_z::T = -0.20
    grid_power::Int = 3
    borrowing_limit_scale::T = 1.0e-7
    asset_max_scale::T = 999.0
    current_account_initial::T = 0.074
    credit_persistence::T = 0.98
    transition_horizon::Int = 200
    egm_tol::T = 1.0e-12
    egm_maxiter::Int = 5000
    dist_tol::T = 1.0e-11
    dist_maxiter::Int = 40000
    root_tol::T = 1.0e-8
    root_maxiter::Int = 80
end

default_parameters(;kwargs...) = validate_params(ModelParams{Float64}(;kwargs...))

function validate_params(p::ModelParams)
    0 < p.alpha < 1 || throw(ArgumentError("alpha must lie in (0,1)"))
    0 < p.delta < 1 || throw(ArgumentError("delta must lie in (0,1)"))
    0 < p.chi < 1 || throw(ArgumentError("chi must lie in (0,1)"))
    p.gamma > 0 || throw(ArgumentError("gamma must be positive"))
    # The implemented c_H Euler equation requires theta = 1/gamma.
    abs(p.theta-one(p.theta)) < 1e-12 && abs(p.gamma-one(p.gamma)) < 1e-12 ||
        throw(ArgumentError("Stage 1 aggregate labor closure requires theta=gamma=1"))
    p.theta_star > 0 || throw(ArgumentError("theta_star must be positive"))
    p.epsilon_f > 1 && p.epsilon_w > 1 ||
        throw(ArgumentError("elasticities must exceed one"))
    p.varphi > 0 || throw(ArgumentError("labor curvature must be positive"))
    abs(p.rho_s) < 1 || throw(ArgumentError("rho_s must be stationary"))
    p.sigma_log_stationary > 0 || throw(ArgumentError("sigma must be positive"))
    p.n_s >= 2 || throw(ArgumentError("n_s must exceed one"))
    p.n_coarse >= 5 && p.n_fine > p.n_coarse ||
        throw(ArgumentError("invalid asset grid sizes"))
    p.grid_min_z < 0 || throw(ArgumentError("negative z-grid portion required"))
    isodd(p.grid_power) || throw(ArgumentError("grid_power must be odd"))
    p.asset_max_scale > 1 || throw(ArgumentError("asset_max_scale too small"))
    p.egm_tol > 0 && p.dist_tol > 0 && p.root_tol > 0 ||
        throw(ArgumentError("tolerances must be positive"))
    return p
end

"""
Reference moments exactly as stated in the publicly available Main.m.
GDP is measured per quarter, while the wealth/NFA inputs are per capita.
"""
struct ReferenceCalibration{T<:AbstractFloat}
    gross_capital::T
    quarterly_gdp::T
    target_nfa::T
    target_wealth::T
    k_y::T
    nfa_y::T
    K::T
    Y::T
    L::T
    r::T
    w::T
    firm_profits::T
    chi_dis::T
    export_shift::T
    ch_budget::T
end

function baseline_reference(p::ModelParams{T}) where {T}
    gdp = T(12.65)/4
    net_wealth = T(26.2)
    nfa = -T(6.5)*T(0.69)
    capital_pc = net_wealth - nfa
    ky = capital_pc/gdp
    K = ky^(inv(1-p.alpha))
    L = one(T)
    Y = K^p.alpha
    markup_inv = (p.epsilon_f-1)/p.epsilon_f
    r = markup_inv*p.alpha*Y/K - p.delta
    w = markup_inv*(1-p.alpha)*Y/L
    profits = Y/p.epsilon_f
    B = (nfa/gdp)*Y
    ch = p.chi*(r*(K+B) + w*L + profits)
    ch > 0 || throw(DomainError(ch, "baseline household consumption"))
    union_markup_inv = (p.epsilon_w-1)/p.epsilon_w
    chi_dis = p.chi*union_markup_inv*w/(ch*L^p.varphi)
    D = Y - p.delta*K - ch
    D > 0 || throw(DomainError(D, "baseline positive foreign demand"))
    return ReferenceCalibration(capital_pc,gdp,nfa,net_wealth,ky,nfa/gdp,
                                K,Y,L,r,w,profits,chi_dis,D,ch)
end
