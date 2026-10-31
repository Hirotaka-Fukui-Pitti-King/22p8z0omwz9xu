#!/usr/bin/env julia
# Independent small-N oracle for the finite-firm switching sweep.
# This script is deliberately O(N^2) and must never be used in simulation.

using Nirei2015

function brute_fixed_points(state, calibration, workspace)
    p = calibration.parameters
    rho = calibration.constants.rho
    n = length(state.lambda)
    boundaries = sort!(unique(vcat(0.0, workspace.switches, Inf)))
    trial_capital = Vector{Float64}(undef, n)
    ks = Float64[]
    xs = Float64[]

    for j in 1:(length(boundaries)-1)
        lower, upper = boundaries[j], boundaries[j+1]
        probe = isfinite(upper) ? lower + (upper - lower) / 2 :
                lower + max(1.0, lower)
        total_investment = 0.0
        for i in 1:n
            d = (1.0 - p.delta) * state.capital[i]
            lambda = state.lambda[i]
            coefficient = threshold_coefficient(
                lambda, state.next_productivity[i], calibration
            )
            newk = if probe <= d / (lambda * coefficient)
                d / lambda
            elseif probe <= d / coefficient
                d
            else
                d * lambda
            end
            trial_capital[i] = newk
            total_investment += newk - d
        end
        K = aggregate_capital(
            state.next_productivity, trial_capital, p.alpha, rho
        )
        if lower < K <= upper
            push!(ks, K)
            push!(xs, total_investment / n)
        end
    end
    return ks, xs
end

function main()
    calibration = run_calibration()
    for seed in (20261101, 20261102, 20261103, 20261104)
        rngs = firm_rngs(seed)
        firms = initialize_firms(calibration, rngs; N=48)
        draw_productivity!(firms, rngs, calibration)
        workspace = initialize_equilibrium_workspace(firms)
        before = copy(firms.capital)
        results = enumerate_fixed_points!(workspace, firms, calibration)
        oracle_ks, oracle_xs = brute_fixed_points(firms, calibration, workspace)

        @assert before == firms.capital
        @assert results.total == length(oracle_ks)
        @assert workspace.n_fixed == length(oracle_ks)
        @assert results.admissible == count(>(0.0), oracle_xs)
        for i in eachindex(oracle_ks)
            @assert isapprox(workspace.fixed_capital[i], oracle_ks[i];
                             atol=1e-10, rtol=1e-11)
            @assert isapprox(workspace.fixed_investment[i], oracle_xs[i];
                             atol=1e-11, rtol=1e-11)
        end
        println("Oracle passed for seed=", seed,
                ", N=48, all fixed points=", length(oracle_ks))
    end
    println("Independent finite-N enumeration oracle passed.")
    return nothing
end

main()
