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

    println("PASS: Strict Section 3.4 rejects nonpositive-only investment roots.")
    println("PASS: Capital fallback selects closest K from all negative-X roots.")
    println("PASS: Positive-X original selection is bit-for-bit unchanged.")
    println("Note: This tests only ES1 selection, not a nonlinear dynamic RE equilibrium.")
end

main()
