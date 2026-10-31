#!/usr/bin/env julia
# Finite-N mechanism audit for the first investment collapse; not a model alternative.
# From Nirei-2015/: julia --project=. scripts/audit_investment_dynamics.jl
using Nirei2015
using TOML

function response_at_K(firms, calibration, K::Float64)
    p = calibration.parameters
    rho = calibration.constants.rho
    n = length(firms.capital)
    power = 0.0
    investment = 0.0
    ninvest = 0
    ninactive = 0
    ndivest = 0
    @inbounds for i in 1:n
        d = (1.0 - p.delta) * firms.capital[i]
        l = firms.lambda[i]
        a = firms.next_productivity[i]
        b = threshold_coefficient(l, a, calibration)
        nextk = if K <= d / (l*b)
            ndivest += 1
            d / l
        elseif K <= d / b
            ninactive += 1
            d
        else
            ninvest += 1
            d*l
        end
        power += a^(rho/p.alpha) * nextk^rho
        investment += nextk - d
    end
    return (
        Phi=(power/n)^(1/rho),
        X=investment/n,
        investing=ninvest,
        inactive=ninactive,
        divesting=ndivest,
    )
end

function gap_diagnostics(state)
    firms, calibration = state.firms, state.calibration
    n = length(firms.capital)
    total = 0.0
    total2 = 0.0
    outsider = 0
    high = -Inf
    low = Inf
    # Current capital and productivity use the same time index.
    @inbounds for i in 1:n
        b = threshold_coefficient(
            firms.lambda[i], firms.productivity[i], calibration
        )
        g = log(firms.capital[i]/(b*state.K)) / log(firms.lambda[i])
        low = min(low, g)
        high = max(high, g)
        outsider += (g < 0.0 || g >= 1.0)
        total += g
        total2 += g^2
    end
    m = total/n
    return Dict(
        "mean_gap"=>m,
        "sd_gap"=>sqrt(max(0.0, total2/n-m^2)),
        "min_gap"=>low, "max_gap"=>high,
        "fraction_gap_outside_0_1"=>outsider/n,
    )
end

function main()
    state = initialize_simulation(; N=110000, base_seed=20261101)
    cal = state.calibration
    p = cal.parameters
    s = cal.steady_state
    coefficient = cal.threshold.A^(cal.constants.rho-1.0)/cal.constants.rho
    n = length(state.firms.capital)

    qsum = 0.0
    violated = 0
    for l in state.firms.lambda
        q = -log1p(-p.delta)/log(l)
        qsum += q
        violated += (q > 1.0)
    end
    output = Dict{String,Any}(
        "configuration"=>Dict(
            "firms"=>n, "seed"=>20261101,
            "continuum_K"=>s.K, "continuum_X"=>s.X,
            "expected_stationary_investment_participation_ignoring_productivity"=>qsum/n,
            "count_lambda_does_not_satisfy_one_sided_condition"=>violated,
            "equation18_investment_multiplier"=>coefficient,
        ),
    )
    periods = Vector{Dict{String,Any}}()
    for t in 0:5
        K_prior = state.K
        expected = compute_expected_equilibrium(
            state.re, cal, K_prior, state.mu_prior
        )
        initial_gap = gap_diagnostics(state)
        draw_productivity!(state.firms, state.rngs, cal)
        enum = enumerate_fixed_points!(state.workspace, state.firms, cal)
        Ke = (1.0-p.delta)*K_prior + coefficient*expected.X_next
        at_expected = response_at_K(state.firms, cal, Ke)
        info = Dict{String,Any}(
            "t"=>t, "K_prior"=>K_prior,
            "expected_K"=>Ke, "expected_X"=>expected.X_next,
            "reaction_at_expected_K"=>at_expected.Phi,
            "reaction_residual_at_expected_K"=>at_expected.Phi-Ke,
            "investment_at_expected_K"=>at_expected.X,
            "investing_at_expected_K"=>at_expected.investing,
            "inactive_at_expected_K"=>at_expected.inactive,
            "divesting_at_expected_K"=>at_expected.divesting,
            "n_fixed"=>enum.total,
            "n_positive_X"=>enum.admissible,
            "min_fixed_K"=>enum.minimum_K, "max_fixed_K"=>enum.maximum_K,
            "min_fixed_X"=>minimum(view(state.workspace.fixed_investment,1:enum.total)),
            "max_fixed_X"=>maximum(view(state.workspace.fixed_investment,1:enum.total)),
            "current_gap"=>initial_gap,
        )
        if enum.admissible == 0
            info["status"] = "no positive investment fixed point; do not change ES1"
            push!(periods, info)
            break
        end
        selection = select_es1(state.workspace,expected.X_next,Ke)
        chosenK = state.workspace.fixed_capital[selection.index]
        chosenX = state.workspace.fixed_investment[selection.index]
        info["status"] = "selected"
        info["K_selected"] = chosenK
        info["X_selected"] = chosenX
        info["epsilon"] = log(chosenX/expected.X_next)
        info["equation18_residual"] = chosenK -
            ((1.0-p.delta)*K_prior + coefficient*chosenX)
        chosen = solve_equilibrium!(
            state.workspace, state.firms, cal, expected.X_next, K_prior
        )
        @assert isapprox(chosen.K_next, chosenK;atol=1e-10,rtol=1e-12)
        state.K = chosen.K_next
        # Existing Issue #5 first-order closure, deliberately unchanged.
        state.mu_prior = s.mu * exp(state.re.h*log(state.K/s.K))
        state.period += 1
        push!(periods, info)
    end
    output["periods"] = periods
    path = normpath(joinpath(@__DIR__,"..","results","initial_collapse_audit.toml"))
    mkpath(dirname(path))
    open(path,"w") do io
        TOML.print(io,output)
    end
    println("Audit saved: ",path)
    println("One-sided-condition violations: ",violated," of ",n)
    for rec in periods
        println("t=",rec["t"]," expected K=",rec["expected_K"],
            " Phi(expected K)-expected K=",rec["reaction_residual_at_expected_K"],
            " selected X=",get(rec,"X_selected","no positive-X fixed point"),
            " fixed points=",rec["n_fixed"])
    end
    return nothing
end

main()
