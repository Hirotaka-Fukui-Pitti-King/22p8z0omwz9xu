# Stage 3 — Source-order repair, 2026-11-03

## Verified from owner Windows logs

- Julia 1.12.6: 34 Julia files parsed, package load passed, Stage 1 2420/2420, Stage 2 37/37, Stage 3 39/39 tests passed.
- Recomputed foreign-debt expansion: scaled GE residual 3.55103e-7.
- **Direct leveraged** surprise-reversal GE did **not converge** through at least iteration 46 and 24,344 residual evaluations. Its maximum scaled residual remained ~0.016 (required 1e-6), concentrated at the first two labor conditions, residual indices 200 and 201.
- Fresh Jacobians and LM fallback did not resolve the stalled direct leveraged formulation. Do not continue the obsolete root solve or represent it as a verified equilibrium.

## Original MATLAB sequence (Main.m lines 164, 342–351; MIT_transition.m 127–154)

1. Compute foreign-denominated credit expansion at k_const_hat=0.
2. Solve **unleveraged** surprise credit reversal at k_const_hat=0.
3. Set k_const_hat=1/16 and **warm-start** the leveraged surprise credit reversal from the converged unleveraged K/L path.
4. Warm-start the fixed-exchange-rate equilibrium from the leveraged flexible solution.

Previous Julia Stage 3 skipped step 2 and attempted step 3 directly from a stationary interpolation. That is a real source-protocol error, not merely a choice of nonlinear solver.

## New committed implementation

- scripts/run_sudden_stop.jl now implements the original two-stage surprise-reversal solution sequence and writes results/stage3_unleveraged_transition.csv **only if the unleveraged equilibrium converges and passes independent diagnostics**.
- scripts/diagnose_stage3_initial.jl performs a one-time residual evaluation for k_const_hat=0 and 1/16 at an identical initial aggregate guess, reporting the FX root, first five labor and capital residuals and initial-market quantities. It does not run any full 399-variable nonlinear solver.
- Original economic equations, borrowing bounds and convergence tolerance remain unchanged.

## Verification limitation

This new patch has been read back from the current branch, but **not syntax-checked or executed in Julia**. Whether the unleveraged 399-equation equilibrium converges is not yet known. The leveraged 399-equation and fixed 599-equation solutions, and Figures 1/4/5, remain unverified. Run scripts/check_syntax.jl and test/runtests.jl on Windows before running scripts/diagnose_stage3_initial.jl. Do not run the expensive full script again unless the one-pass diagnostic confirms the initial conditions are feasible.
