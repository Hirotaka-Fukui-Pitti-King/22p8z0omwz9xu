using Nirei2015

"""
Single, deliberately exogenous-K one-period test of the stationary-gap
assumption in Nirei (2015), Section 2.6 and Appendix C.

Keep both firm policies and productivity draws exactly as in the original
three-choice code. This does NOT solve a finite-N GE fixed point or ES1.
Run locally: julia --project=. scripts/check_stationary_one_step.jl 110000 20261101
"""
function check_stationary_one_step(N::Int, seed::Int)
    N > 0 || throw(ArgumentError("N must be positive"))
    seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
    cal = run_calibration()
    rngs = firm_rngs(seed, 1)
    firms = initialize_firms(cal, rngs; N=N)
    p = cal.parameters
    Kbar = cal.steady_state.K
    K0 = aggregate_capital(firms, cal)
    draw_productivity!(firms, rngs, cal)

    k1_no_shock = Vector{Float64}(undef, N)
    k1_shock = Vector{Float64}(undef, N)
    outside_no_shock = 0
    outside_shock = 0
    below_shock = 0
    above_shock = 0
    lambda_support_violation = 0
    n_investing = 0
    n_divesting = 0
    gap_identity_max_error = 0.0
    investment_no_shock = 0.0
    investment_shock = 0.0

    @inbounds for i in 1:N
        lambda_i = firms.lambda[i]
        oldk = firms.capital[i]
        dcap = (1.0 - p.delta) * oldk
        b0 = threshold_coefficient(lambda_i, firms.productivity[i], cal)
        b1 = threshold_coefficient(lambda_i, firms.next_productivity[i], cal)
        lambda_support_violation += lambda_i * (1.0 - p.delta) <= 1.0

        # Thresholds are the original published {−1,0,+1} policy, with
        # exogenous target aggregate capital fixed at the continuum Kbar.
        action0 = dcap < b0 * Kbar ? 1 :
            (dcap >= lambda_i * b0 * Kbar ? -1 : 0)
        action1 = dcap < b1 * Kbar ? 1 :
            (dcap >= lambda_i * b1 * Kbar ? -1 : 0)

        k_no_shock = dcap * lambda_i^action0
        k_shock = dcap * lambda_i^action1
        k1_no_shock[i] = k_no_shock
        k1_shock[i] = k_shock
        investment_no_shock += (k_no_shock - dcap)
        investment_shock += (k_shock - dcap)

        s_no_shock = log(k_no_shock / (b0 * Kbar)) / log(lambda_i)
        s_shock = log(k_shock / (b1 * Kbar)) / log(lambda_i)
        outside_no_shock += (s_no_shock < 0.0 || s_no_shock >= 1.0)
        outside_shock += (s_shock < 0.0 || s_shock >= 1.0)
        below_shock += (s_shock < 0.0)
        above_shock += (s_shock >= 1.0)
        n_investing += (action1 == 1)
        n_divesting += (action1 == -1)

        # Independently verify the exact gap-transition identity, including
        # the productivity-induced threshold displacement.
        predicted_gap = firms.initial_gap[i] +
            (log1p(-p.delta) + log(b0 / b1)) / log(lambda_i) + action1
        gap_identity_max_error = max(gap_identity_max_error,
                                     abs(predicted_gap - s_shock))
    end

    K1_no_shock = aggregate_capital(
        firms.productivity, k1_no_shock, p.alpha, cal.constants.rho
    )
    K1_shock = aggregate_capital(
        firms.next_productivity, k1_shock, p.alpha, cal.constants.rho
    )

    println("Nirei (2015) one-step *exogenous-K* stationarity experiment")
    println("N=", N, " seed=", seed, " Kbar=", Kbar)
    println("Initial finite-N CES K0/Kbar=", K0 / Kbar)
    println("Lambda support violation fraction=", lambda_support_violation / N)
    println("No-productivity-shock gap outside fraction=", outside_no_shock / N)
    println("IID productivity-shock gap outside fraction=", outside_shock / N)
    println("  below zero=", below_shock / N, " above one=", above_shock / N)
    println("IID shock investing share=", n_investing / N,
            " divesting share=", n_divesting / N)
    println("No-shock K1/K0-1=", K1_no_shock / K0 - 1.0,
            " investment=", investment_no_shock / N)
    println("Shock K1/K0-1=", K1_shock / K0 - 1.0,
            " investment=", investment_shock / N)
    println("Gap transition max algebraic residual=", gap_identity_max_error)
    gap_identity_max_error <= 1e-9 ||
        error("Exact gap transition algebra failed: $gap_identity_max_error")
    println("This DOES NOT establish a full GE equilibrium or Table 2 replication.")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 110_000
    seed = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20261101
    check_stationary_one_step(n, seed)
end
