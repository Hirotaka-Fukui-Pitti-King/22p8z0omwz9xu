# FMR-2020 Stage 3 — Author source audit and implementation status

## Canonical source and stage sequencing

Read the original Main.m, MIT_transition.m, Transition/solve_transition_bis.m, Transition/solve_transition_fixed.m, Transition/backsolve_egm.m, and Figures.m before implementing. The original Main.m sets KK=2 so `fixed=false, foreign=true` in MIT_transition.m. In this experiment all households receive the original first-period debt revaluation even if aggregate B_initial=0. When k_const_hat=0 and B_initial=0, the code sets leverage=1. Positive-wealth households then have no gross foreign bonds; negative-wealth households can retain a tiny net foreign liability. Consequently the previous Stage 2 baseline was not **literally** identical and must be recalculated with the correct source convention. The Stage 3 script does this.

The original source switches `params.k_const_hat=1/16` before solving the sudden stop and the fixed-regime counterfactual. The reversal `shock_start_t=41` inherits the original period-41 capital and wealth distribution, but the foreign bond balance B and capital price q from period 40. The `terms_initial` is also the pre-shock period-40 value, not the post-shock price. A source-code-level output table is required to inspect the implied wealth before and after revaluation.

## Household balance sheets

For `a` on the full fine grid, compute

`K_h=max(lambda*a + k_const,0)`,
`B_h=a-K_h`,
`a_new=K_h+B_h/(terms_new/terms_pre)`.

`k_const=k_const_hat*K_initial`. With nonzero inherited foreign debt `lambda=(1-k_const_hat)*K_initial*q_initial/rich_net_wealth`, where `rich_net_wealth=sum_{a>0,s}mu(a,s)*a`. With B_initial=0, the author code overrides lambda=1. Full fine-grid negative states are retained to accommodate foreign liability revaluation, rather than projecting them onto the stationary borrowing limit. Probability is split between neighboring asset knots. If any positive-mass revalued wealth falls outside this full grid or has negative/nonfinite feasible consumption, the code raises a diagnostic error, not clipping.

The flexible first-period trade-price scalar condition uses aggregate opening assets `q_initial*K_initial + B_initial*terms_pre/terms_1`. This matches solve_transition_bis.m for foreign=true. Other dates use standard Stage 2 wealth. The first-period household distribution is revalued after terms_1 is calculated, so every proposed aggregate path corresponds to its own cross-sectional balance sheets.

## 399-variable flexible and 599-variable fixed

Flexible: unknown `[K[2:T];L[1:T]]`, equation count `(T-1)+T=399`. The original term price is solved from the domestic goods budget; backward EGM, forward wealth distribution, asset and labor clearing are recomputed per residual.

Fixed: unknown `[K[2:T];L[1:T];terms[1:T]]`, equation count `(T-1)+T+T=599`. Relative domestic price inflation is `pi_1=terms_1/terms_pre, pi_t=terms_t/terms_{t-1}`. Backward real marginal costs obey

`mc_T=(eps_f-1)/eps_f +(theta_nr/eps_f)*(pi_T-1)*pi_T`

`mc_t=(eps_f-1)/eps_f +(theta_nr/eps_f)*(pi_t-1)*pi_t - (theta_nr/eps_f)*(pi_{t+1}-1)*pi_{t+1}*Y_{t+1}/[(1+r_{t+1})Y_t]`.

`w_t=mc_t(1-alpha)Y_t/L_t`, `rK_t=mc_t*alpha*Y_t/K_t`, and `r_t=[rK_t+(1-delta)q_t]/q_{t-1}-1`.

Firm profits are `Y_t-rK_t*K_t-w_t*L_t-(theta_nr/2)(pi_t-1)^2` plus investment-firm profits. Fixed-regime goods clearing includes the negative resource cost `(theta_nr/2)(pi_t-1)^2`, real capital adjustment costs and foreign demand `D*terms_t^(-theta_star)`. The final period follows the author's penultimate investment convention and duplicates consumption from T-1 for the labor and goods terminal residuals.

### Original source discrepancies, made explicit

- **Fixed `solve_transition_fixed.m`, lines 64–67**: the investment adjustment-cost expression in the interior goods equation uses stale `t` after an earlier loop, and its denominators differ from the correctly indexed `i_agg(t)/k_agg(t)` expression. Moreover it contains `- -params.theta_nr/2` instead of subtracting the positive price adjustment resource cost. Stage 3 Julia implements the internally consistent period-t capital adjustment costs and *negative* nominal rigidity resource costs. This is an economic correction; it is not a byte-for-byte numerical replication of that buggy source line.
- **Figure 5, line 351**: the original code computes inverse terms of trade but subtracts `log(tot0)`, rather than `log(1/tot0)`. Stage 3 uses the theoretically consistent change in inverse terms of trade, documenting the source difference.
- **Stage 2 author labor condition**: transition equation omits `1/varphi` power present in stationary union FOC. This source inconsistency is retained, independently logged, and cannot be presented as a resolved economic identification.

### Source files and scope

New: `src/foreign_debt_revaluation.jl`, `src/sudden_stop_equilibrium.jl`, `src/transition_fixed_exchange.jl`, `src/stage3_figures.jl`, `scripts/run_sudden_stop.jl`, `scripts/validate_sudden_stop.jl`, `scripts/plot_figures_1_4_5.jl`, `test/test_stage3.jl`. Modified `src/FMR2020.jl` exports/includes, `test/runtests.jl` includes Stage3 tests, and docs/README/policy.

All Julia source is freshly authored, and **no Julia execution is performed by the assistant**. The owner runs Julia 1.12.6 on Windows. At time of commit: syntax unverified, Stage 3 unit tests unverified, smoke test unverified, 399- and 599-dimensional GE numerical convergence unverified, SVG generation unverified and numerical correspondence with author figures unverified.

Previous Stage 2 owner-PC verification remains true only for the Stage 2 source, and cannot be reused as proof that Stage 3 works. A completed source implementation is not a verified replication. If the owner's first smoke test fails on distribution support or the fixed Rotemberg block, the user-PC stack trace controls the precise subsequent fix; there is no warrant for loosening bounds or claiming success.

## 2026-11-03 owner-PC first-run failure: first-period foreign-debt trade root

Owner ran Julia 1.12.6 on Windows: 34/34 Julia files parsed, 2,420 Stage 1, 37 Stage 2 and 15 Stage 3 tests passed; no-revaluation and fixed-price-block smoke passed (599 unknowns). Recalculated foreign-debt expansion converged immediately from the Stage 2 solution (max scaled GE residual 3.55103e-7, max goods level 5.37672e-7). The 399-variable surprise-reversal solver failed **before its first nonlinear iteration** in stage3_positive_trade_root: the implementation required opposite signs at x=0.001*xprev and x=1000*xprev.

On the actual failed initial aggregate guess, the trade-price residual f(x)=C+E/x-D/x^theta_star has approximately C=-2.77822301392, E=3.82570975404, D=1.01494635970, theta_star=3. It has a maximum at x_peak=0.89212593177 with f_peak=+0.08064822670, and two genuine positive roots, x_low=0.78658958396 and x_high=1.03791664536. Both far endpoints have negative f, so the previous code made a mathematically false no-root inference.

The original source solve_transition_bis.m uses MATLAB fsolve with an initial 0.8 in the final price iteration; the observed initial-guesstimate 0.8 lies near x_low and a Newton step leads to the lower root. Correction: compute the unique positive maximum, construct a strict sign-changing bracket only on (0,x_peak), and solve the lower-price branch with the same original budget/resource equation; explicitly throw when no maximum above zero exists. The independent residual and household/distribution constraints stay unchanged. A unit test and short smoke test now cover the exact observed two-positive-root failure. This is a direct defect repair, not a new numerical economic model.

**The source patch is committed but not yet run on the owner's PC.** The macro GE and Figures 1/4/5 have NOT been solved/generated and remain unverified.

## 2026-11-03: blocked module load from misplaced consecutive Julia docstrings

Owner Windows Julia 1.12.6: 34/34 files printed PARSE OK on the previous commit, but `using FMR2020` failed during precompilation with `LoadError: cannot document the following expression` at `src/sudden_stop_equilibrium.jl:2`. The previous change accidentally placed TWO consecutive triple-quoted docstrings at the top of that file, leaving the first documentation expression with no valid function/struct to document. All three subsequent commands (unit tests, Stage 3 smoke and full nonlinear run) consequently failed at package loading; **none of their economics were executed in that attempt**.

Fix (commit `ef781a05905ffc77931076340cbfd442bb2fe1ca`): convert the first module-level descriptive docstring to comments; the second triple-quoted string now immediately precedes the intended `foreign_debt_trade_root` function. Scanned all 19 `src/*.jl` files for adjacent triple-quote blocks; only this one occurrence was found. Modified `scripts/check_syntax.jl` to **actually load FMR2020 after Meta.parseall**, so future malformed docstring attachments fail the first check instead of displaying false success. Existing parse-only output alone was inadequate; no change to the model, root choice, constraints or numerical tolerances.

**This corrected commit has NOT been tested on the owner's Windows PC.** Until new PC logs, Julia syntax+module-load, full unit tests and macro GE convergence are unverified for the current commit. Do not call Stage 3 completed merely because source was committed.


## 2026-11-03: Positive-mass household budget infeasibility after foreign-liability revaluation

Owner's latest Windows 11, Julia 1.12.6 log: **34/34 parse+package load, 2,420 Stage 1 + 37 Stage 2 + 22 Stage 3 tests passed**, smoke passed (no-revaluation mass test, double-root regression, fixed price block). The foreign-denominated expansion already passes its GE aggregate tests (scaled residual 3.5510e-7). The sudden-stop flexible transition fails before the first nonlinear GE iteration:

`DomainError with (1, -0.6604250433394604, 1, NaN, NaN): positive-mass household infeasible after FX shock`.

This occurs because the previous lower-branch trade root `x_low≈0.78659` revalues gross foreign liabilities sharply upward and sends positive household mass into wealth-grid nodes where the budget cannot fund positive consumption with `a' >= bmin`. That is an **economically infeasible household state**, not a harmless off-support policy node. The original MATLAB backsolve_egm.m lines 68–71 uses `max(c_constrained,1e-5)`, and lines 100–107 use constrained savings `bmin`; for such a state, blindly copying the floor would violate the budget constraint. We deliberately did NOT insert an artificial consumption floor or drop distribution mass.

The trade-price residual has **two** positive roots at this initial aggregate path. The lower root cannot automatically be equated with the correctly selected equilibrium merely because MATLAB `fsolve` started at 0.8. The modified algorithm computes both positive roots and uses the original revaluation mapping plus *the actual post-revaluation mass-splitting grid* to test a necessary household feasibility condition at every occupied (asset,productivity) state:

`(1+r_1)*a_revalued_grid + w_1 L_1 s + profits_1 - bmin > 0`.

It selects the lower root only if that branch satisfies the household budget condition; otherwise it uses the upper root if that branch is feasible. An out-of-grid positive-mass revaluation rejects just that root, not every possible root. If neither branch is feasible, it fails with both roots and the minimum available household resource diagnostic. The later backward EGM and forward distribution still enforce positive consumption, the borrowing lower bound, probability conservation and independently checked household budgets. **No tolerances, structural equations, shock timing or wealth values are altered.** Root selection is now conditional on household feasibility and thus **is not asserted to be numerically identical** to MATLAB's unguarded `fsolve(aggr,0.8)` with consumption floors.

Important: this code patch has been statically reviewed only and has **not** yet passed owner-PC Julia syntax/runtime or shown successful 399/599 GE convergence. Full Stage 3 results remain unverified.
