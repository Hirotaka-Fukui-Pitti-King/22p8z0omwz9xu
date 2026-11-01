
using FMR2020

function main()
    println("Stage 3: validate source shock indices, revaluation and 599-variable price block")
    p=default_parameters(n_coarse=61,n_fine=241)
    st=solve_three_steady_states(p)
    ctx=stage3_expansion_context(st)
    n=length(ctx.economy.B)
    reversal=credit_reversal_path(ctx.economy.B,st.calibrated.B)
    @assert length(reversal)==200
    @assert reversal[1]==ctx.economy.B[40]
    @assert abs(reversal[end]-st.calibrated.B)<1e-11
    # Revaluation with unchanged exchange rate must leave wealth unchanged
    # except for floating-point mass interpolation errors.
    before=ctx.mu_pre
    after=similar(before)
    grid=ctx.economy.grids.fine
    diagnosis=revalue_foreign_debt!(after,before,grid,st.initial.K,1.0,
                                   0.0,st.initial.terms,st.initial.terms,0.0)
    @assert abs(diagnosis.mass_after-1)<1e-10
    @assert abs(diagnosis.net_wealth_after-diagnosis.net_wealth_before)<1e-8
    println("No-revaluation test: ",diagnosis)
    # Fixed-rate nominal Phillips curve in a stationary path:
    # no inflation implies constant markup marginal cost.
    inp=ctx.economy
    x=initial_transition_guess(inp)
    fixedguess=vcat(x,fill(st.initial.terms,n))
    workspace=TransitionWorkspace(inp)
    FMR2020.stage3_fixed_prices!(workspace,ctx,fixedguess)
    @assert all(isfinite,workspace.rate)
    @assert minimum(workspace.terms)>0
    @assert length(fixedguess)==3*n-1
    println("Fixed-rate price block PASS, dimension=",length(fixedguess))
    println("STAGE3 SMOKE PASSED; full GE remains unverified until separate solve.")
end

main()
