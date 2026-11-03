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

## 2026-11-03: Owner one-pass residual results and low-cost staged solve

Owner Windows Julia 1.12.6: 35 files parsed, module load succeeded, 2420+37+39=2496 unit tests passed. Stage 3 corrected foreign-denominated credit expansion verified at maximum scaled GE residual 3.551030838266911e-7.

The new one-pass `scripts/diagnose_stage3_initial.jl` completed successfully. Both cases have strictly positive minimum resources in each occupied revalued household state, and probability mass is conserved. These tests verify initial feasible household budgets, NOT full equilibrium.

| Statistic at identical original-source initial K/L guess | No leverage | Leverage 1/16 |
|---|---:|---:|
| Selected initial terms of trade | 0.7865895866917804 | 1.037916643090368 |
| Min resources across occupied post-FX nodes | 0.4708108248444556 | 0.48493043012136894 |
| Full GE residual 2-norm | 0.4992912129602136 | 0.14595616021482744 |
| Max absolute scaled asset residual | 0.05651597863384016 | 0.009992184169123416 |
| Max absolute scaled labor residual | 0.0420581405452537 | 0.021003904982795363 |
| Period-1 labor residual | 0.0420581405452537 | 0.021003904982795363 |

**Interpretation:** The no-leverage starting vector is further from GE than the leveraged vector in norm. Neither norm establishes existence or nonexistence of a fully feasible GE. The original MATLAB Main.m calls fsolve on no-leverage first (MIT_transition.m), then leverages warm-start from its solution. Julia has been changed to do the same, but this new solve is unverified.

**Process safety update:** `scripts/run_sudden_stop.jl --unleveraged-only` now computes only the unleveraged 399-variable reversal, independently verifies all GE and household diagnostics, saves `results/stage3_unleveraged_transition.csv` and exits. The full run can reuse that saved K/L vector as a warm start, rechecking every GE residual before accepting it. `stage3_nonlinear_solve` now enforces at most 5000 exact residual evaluations per nonlinear solve, matching original MATLAB `MIT_transition.m` `MaxFunEvals=5000`. The evaluation cap is a time/cost guard, not a relaxation of equilibrium conditions. On failure it throws and DOES NOT produce plots as if a solution existed.

The new source changes are not yet Julia syntax-tested, loaded, unit-tested or numerically run on the owner's PC. The earlier 2496 tests belong only to the previous commit. Do not claim completion of any Stage 3 nonlinear GE.


## 2026-11-03: unleveraged full GE hit evaluation limit (owner Windows Julia)

The owner ran Julia 1.12.6 on the prior source: all **35** Julia files parsed and FMR2020 loaded; **2,420 + 37 + 39 = 2,496** unit tests passed. The original-source-order unleveraged credit reversal (initial positive household-resource min 0.4708108) proceeded but **did NOT solve the 399 equations within the 5000 evaluation cap**:

- iteration 0: max scaled residual **0.05651597863384016** and residual 2-norm **0.4992912129602136**;
- iteration 11: max scaled residual **0.00897660738511874** and 2-norm **0.03254463187358031**;
- limit hit while rebuilding a 399-column Jacobian at iteration 12. This is evidence of substantial initial improvement but not convergence and not proof of no solution.
- The first few accepted steps lowered the 2-norm by 58.7%, 44.0%, and 7.65%. Repeated expensive numerical Jacobian recalculations on small step sizes consumed the function budget.

**Numerical work-cost correction**: drop the policy that refreshes all 399 Jacobian columns simply because an accepted line-search damping was <1/8. Instead, refresh when two consecutive improvements are <2%, on the scheduled refresh interval, or when an existing Broyden direction fails. Also avoid a second redundant refresh in the same iteration. Keep the same economic residuals and strict 1e-6 acceptance check.

**Checkpointing**: `stage3_write_checkpoint` writes only accepted but **UNCONVERGED** 399-vector states to a distinct file `results/stage3_unleveraged_UNCONVERGED.checkpoint` with a signature, dimension, iteration, maximum residual and values. It is NEVER treated as `stage3_unleveraged_transition.csv` or accepted output. Subsequent `--unleveraged-only` runs start from that checkpoint if no validated CSV exists, recomputing all current residuals; on verified convergence the checkpoint is removed and the validated CSV is written. Atomic temporary-file replacement avoids partial writes. A small unit test checks roundtrip, dimensions and finite values. The 5000 residual evaluations per run are retained.

**Verification:** This new patch has been statically reviewed only, with no duplicate consecutive docstrings and correct LM function docstring attachment. The owner has not run syntax/load, test or 399-variable GE on the latest commit. Neither the 399-variable flexible sudden stop nor the 599-variable fixed sudden stop is verified.
