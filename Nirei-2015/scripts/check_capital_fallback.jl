#!/usr/bin/env julia
# Select-only regression for original Section 3.4 and opt-in Section 2.3
# capital-based fallback. No GE simulation, shock rejection or new policies.
# From Nirei-2015/: julia --project=. scripts/check_capital_fallback.jl

using Nirei2015

function main()
    workspace = initialize_equilibrium_workspace(2)
    workspace.n_fixed = 3
    workspace.fixed_capital[1:3] .= [7.90, 8.02, 8.08]
    workspace.fixed_investment[1:3] .= [-0.20, -0.08, -0.02]

    rejected = false
    try
        select_es1(workspace, 0.17, 8.07)
    catch err
        rejected = err isa DomainError
    end
    @assert rejected "Original Section 3.4 must reject all-negative X"

    fallback = select_es1(
        workspace, 0.17, 8.07; selection_mode=:capital_fallback
    )
    @assert fallback.index == 3
    @assert fallback.admissible == 0
    @assert fallback.nonpositive == 3
    @assert isapprox(fallback.criterion, abs(log(8.08 / 8.07)); atol=1e-14)

    # The explicit completion must not alter Section 3.4 if an X>0 root
    # exists; comparison of indices AND objectives catches this.
    workspace.fixed_investment[2] = 0.15
    original = select_es1(workspace, 0.17, 8.07)
    completed = select_es1(
        workspace, 0.17, 8.07; selection_mode=:capital_fallback
    )
    @assert original.index == completed.index == 2
    @assert original.criterion == completed.criterion
    @assert original.admissible == completed.admissible == 1

    # All-root capital criterion must select the closest capital root
    # even when a positive-X root exists elsewhere.
    capital_all = select_es1(
        workspace, 0.17, 8.07; selection_mode=:capital_all
    )
    @assert capital_all.index == 3
    @assert capital_all.criterion == abs(log(8.08 / 8.07))
    @assert capital_all.admissible == 1
    @assert capital_all.nonpositive == 2
    @assert capital_all.index != original.index

    # With only negative investment, all-root mode still selects a
    # genuine exact capital fixed point; it must never reject for X<0.
    workspace.fixed_investment[2] = -0.08
    all_negative = select_es1(
        workspace, 0.17, 8.07; selection_mode=:capital_all
    )
    @assert all_negative.index == 3
    @assert all_negative.admissible == 0

    # A zero-investment closest root is also a legitimate capital root
    # for the capital selection, though log(X/Xe) is undefined.
    workspace.fixed_investment[3] = 0.0
    zero_root = select_es1(
        workspace, 0.17, 8.07; selection_mode=:capital_all
    )
    @assert zero_root.index == 3
    @assert zero_root.nonpositive == 3

    println("PASS: Strict Section 3.4 rejects nonpositive-only investment roots.")
    println("PASS: Capital fallback selects closest K when X>0 is unavailable.")
    println("PASS: Positive-X original selection is unchanged in fallback mode.")
    println("PASS: Capital_all selects all-root closest K even when a positive-X root exists.")
    println("PASS: Capital_all admits nonpositive and exactly zero investment.")
    println("Note: These tests exercise selection only, not a nonlinear dynamic RE equilibrium.")
end

main()
