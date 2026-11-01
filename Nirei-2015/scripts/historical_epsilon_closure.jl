#!/usr/bin/env julia
# Historical author epsilon/continuum-capital protocol — NOT exact 2015 finite-N GE.
# The firms still use the published 2015 {-1,0,+1} policy, Table 1 sigma=1.5,
# unchanged shocks and exact firm CES reaggregation. In contrast to published
# Eq. (13) as a FIXED POINT, the candidate K is determined from continuum Eq. (18)
# Khat = (1-delta)*Kprior + A^(rho-1)/rho * X(Khat).
# The discrepancy between selected Khat and exact CES of selected firms is
# measured and NOT silently treated as satisfying the published equilibrium.
# From Nirei-2015: julia --project=. scripts/historical_epsilon_closure.jl [T] [N] [seed]

using Nirei2015
using TOML
using Statistics: std

function select_continuum_epsilon!(state::SimulationState, expected_X::Float64)
    firm = state.firms
    cal = state.calibration
    p = cal.parameters
    rho = cal.constants.rho
    n = length(firm.capital)
    firm.next_ready || error("draw next productivity before continuum selection")
    invn = inv(Float64(n))
    capital_coefficient = cal.threshold.A^(rho-1.0) / rho
    baseline_capital = (1.0-p.delta)*state.K
    switches = Vector{Float64}(undef, 2n)
    depreciated = Vector{Float64}(undef, n)
    # Eq. (14): signed investment relative to depreciated existing assets.
    investment = 0.0
    @inbounds for i in 1:n
        d = (1.0-p.delta)*firm.capital[i]
        lambda = firm.lambda[i]
        b = threshold_coefficient(lambda, firm.next_productivity[i], cal)
        depreciated[i] = d
        switches[2i-1] = d/(lambda*b)  # divest -> inactive
        switches[2i] = d/b             # inactive -> invest
        investment += d*(inv(lambda)-1.0)*invn
    end
    order = sortperm(switches)
    left = 0.0
    j = 1
    # The original published Section 3.4 chooses the root nearest Xe in log X.
    best_K = NaN
    best_X = NaN
    best_score = Inf
    n_candidates = 0
    n_positive = 0
    function record_candidate(right::Float64, current_X::Float64)
        candidate = baseline_capital + capital_coefficient*current_X
        if isfinite(candidate) && left < candidate <= right
            return candidate
        end
        return NaN
    end
    while j <= 2n
        boundary = switches[order[j]]
        candidate = record_candidate(boundary, investment)
        if isfinite(candidate)
            n_candidates += 1
            if investment > 0.0
                n_positive += 1
                score = abs(log(investment / expected_X))
                if score < best_score
                    best_K = candidate
                    best_X = investment
                    best_score = score
                end
            end
        end
        # Breakpoint equality uses the PRE-switch policy, as in original code.
        while j <= 2n && switches[order[j]] == boundary
            event = order[j]
            i = (event+1) >>> 1
            lambda = firm.lambda[i]
            d = depreciated[i]
            if isodd(event)
                investment += d*(1.0-inv(lambda))*invn
            else
                investment += d*(lambda-1.0)*invn
            end
            j += 1
        end
        left = boundary
    end
    candidate = record_candidate(Inf, investment)
    if isfinite(candidate)
        n_candidates += 1
        if investment > 0.0
            n_positive += 1
            score = abs(log(investment / expected_X))
            if score < best_score
                best_K = candidate
                best_X = investment
                best_score = score
            end
        end
    end
    isfinite(best_K) && best_X > 0.0 ||
        error("No positive-X continuum Eq18 approximate root (roots=$n_candidates, positive=$n_positive)")
    next_capital = Vector{Float64}(undef, n)
    x_sum = 0.0
    investing = 0
    inactive = 0
    divesting = 0
    @inbounds for i in 1:n
        d = depreciated[i]
        if best_K <= switches[2i-1]
            newk = d/firm.lambda[i]
            divesting += 1
        elseif best_K <= switches[2i]
            newk = d
            inactive += 1
        else
            newk = d*firm.lambda[i]
            investing += 1
        end
        next_capital[i] = newk
        x_sum += newk-d
    end
    X_direct = x_sum*invn
    abs(X_direct-best_X) < 1e-9*max(1.0,abs(best_X)) ||
        error("direct individual investment fails approximate Eq18 policy: $(X_direct-best_X)")
    K_direct = aggregate_capital(firm.next_productivity, next_capital, p.alpha, rho)
    # This is the genuine 2015 Eq13, usually DIFFERENT from Eq18 candidate.
    gap = K_direct-best_K
    firm.capital = next_capital
    advance_productivity!(firm)
    return (Khat=best_K, Kfirm=K_direct, X=best_X, K_gap=gap,
            n_roots=n_candidates, n_positive=n_positive,
            n_investing=investing, n_divesting=divesting,
            n_inactive=inactive, investment_residual=X_direct-best_X,
            epsilon=log(best_X/expected_X))
end

function simulate_approx_period!(state::SimulationState)
    ss = state.calibration.steady_state
    period = state.period
    kp = state.K
    mu_prior = state.mu_prior
    expected = compute_expected_equilibrium(state.re,state.calibration,kp,mu_prior)
    draw_productivity!(state.firms,state.rngs,state.calibration)
    eq = select_continuum_epsilon!(state,expected.X_next)
    # Carry actual CES aggregate forward, not the Eq18 approximation.
    logK_prior = log(kp/ss.K)
    logK_next = log(eq.Kfirm/ss.K)
    logmu = state.re.h*logK_next
    next_mu = ss.mu*exp(logmu)
    macro_vars = Nirei2015._realized_macroeconomy(
        state.calibration,state.re,logK_prior,log(eq.X/ss.X),
        logmu,log(mu_prior/ss.mu))
    level_resource_gap = (macro_vars.Y-macro_vars.C-eq.X)/macro_vars.Y
    state.K = eq.Kfirm
    state.mu_prior = next_mu
    state.period += 1
    return (t=period,Kfirm=eq.Kfirm,Khat=eq.Khat,K_gap=eq.K_gap,
            X=eq.X,X_expected=expected.X_next,epsilon=eq.epsilon,
            n_roots=eq.n_roots,n_positive=eq.n_positive,
            n_investing=eq.n_investing,n_inactive=eq.n_inactive,
            n_divesting=eq.n_divesting,
            direct_X_residual=eq.investment_residual,
            Y=macro_vars.Y,C=macro_vars.C,L=macro_vars.L,
            level_resource_gap=level_resource_gap)
end

function main()
    T = length(ARGS)>=1 ? parse(Int,ARGS[1]) : 20
    N = length(ARGS)>=2 ? parse(Int,ARGS[2]) : 110000
    seed = length(ARGS)>=3 ? parse(Int,ARGS[3]) : 20261101
    T>0 && N>0 || error("T and N must be positive")
    state = initialize_simulation(N=N,base_seed=seed)
    results = NamedTuple[]
    failure = "none"
    for _ in 1:T
        try
            rec = simulate_approx_period!(state)
            push!(results,rec)
            if rec.t<=5 || rec.t==T-1
                println("t=$(rec.t) Khat=$(rec.Khat) K_CES=$(rec.Kfirm) CES_gap=$(rec.K_gap) X=$(rec.X) Xe=$(rec.X_expected) eps=$(rec.epsilon) approximate_roots=$(rec.n_roots)")
            end
        catch ex
            failure=sprint(showerror,ex)
            println(stderr,"Approximate continuum closure stopped t=$(state.period): ",failure)
            break
        end
    end
    dest=normpath(joinpath(@__DIR__,"..","results"))
    mkpath(dest)
    prefix="historical_epsilon_T$(T)_N$(N)_seed$(seed)"
    csv=joinpath(dest,prefix*".csv")
    fields=(:t,:Kfirm,:Khat,:K_gap,:X,:X_expected,:epsilon,
            :n_roots,:n_positive,:n_investing,:n_inactive,:n_divesting,
            :direct_X_residual,:Y,:C,:L,:level_resource_gap)
    open(csv,"w") do io
        println(io,join(string.(fields),","))
        for r in results
            println(io,join((string(getfield(r,k)) for k in fields),","))
        end
    end
    summary=Dict(
        "description"=>"Diagnostic continuum Eq18 epsilon-equilibrium, original 2015 3 choices and sigma=1.5. NOT exact finite-N Eq13 GE.",
        "T"=>T,"N"=>N,"seed"=>seed,
        "completed_periods"=>length(results),
        "completed"=>length(results)==T,"failure"=>failure,
        "max_abs_CES_mismatch"=>isempty(results) ? 0.0 : maximum(r->abs(r.K_gap),results),
        "max_abs_X_reaggregation_residual"=>isempty(results) ? 0.0 : maximum(r->abs(r.direct_X_residual),results),
        "max_abs_level_resource_gap"=>isempty(results) ? 0.0 : maximum(r->abs(r.level_resource_gap),results))
    if length(results)>100
        kept=results[101:end]
        summary["postburnin"]=Dict(
            "Y_sd_percent"=>100*std(log(r.Y) for r in kept),
            "C_sd_percent"=>100*std(log(r.C) for r in kept),
            "X_sd_percent"=>100*std(log(r.X) for r in kept),
            "K_sd_percent"=>100*std(log(r.Kfirm) for r in kept))
    end
    path=joinpath(dest,prefix*".toml")
    open(path,"w") do io
        TOML.print(io,summary)
    end
    println("RESULT: ",summary)
    println("Saved: ",csv," and ",path)
    length(results)==T || error("approximate continuum trial stopped before T at t=$(state.period)")
end

main()
