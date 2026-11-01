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
