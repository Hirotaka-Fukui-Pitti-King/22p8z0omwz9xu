# FMR (2020) source-to-Julia algorithm audit — 2026-11-03

## Reproduction target and observed failure

Original: [kurtmitman/HANKSOME](https://github.com/kurtmitman/HANKSOME).
Julia: [22p8z0omwz9xu / FMR-2020](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/tree/research/fmr-2020-stage3-sudden-stop-exchange-regimes/FMR-2020).
The latest owner-PC run (Windows Julia 1.12.6) parsed and loaded 35 Julia files, and 2420+37+45=2502 tests passed. Credit expansion has max scaled GE residual 3.551030838266911e-7. The **unleveraged sudden-stop 399-dimensional problem did not converge**: at iteration 16, maximum scaled residual 0.0051466392021841875 (asset residual index 1), norm 0.03289858309733816, with 3499 function evaluations. Leveraged flexible and fixed-price equilibria, and Figures 1/4/5, remain unverified.

Do not interpret a small fixed-point update, a stable approximate Jacobian, 2502 passed component tests, or a 5000-call budget as numerical reproduction of the source MATLAB solution.

## Direct line-by-line functional map (verified by reading exact source, NOT by cross-language runs)

| Original MATLAB (file:lines) | Julia (file, functions) | Operation | Status and reproduction consequence |
|---|---|---|---|
| Main.m:95-160 | grids.jl:make_asset_grids | Cubic grid on z in [-0.2,1] with 193/769 nodes | **Different effective grids**: Julia stores full 193/769, but only uses 161 active coarse states for EGM; original operates on full 193 coarse. It matters at least for interpolation nodes, terminal policy dimension and near-constraint behavior. |
| Main.m:125-161 | grids.jl, sudden_stop_equilibrium.jl:stage3_household_backward! | Original Index_b_min = 33, fine = 129, while index 1 is strictly below bmin | **Different EGM cutoff**: Julia formerly tests first active node, whereas source uses endogenous asset at index 33 of full grid. The active policy grids make these economically related, but they are not the same discretization. |
| Steady_State_Ayagari/solve_R.m:18-75 | household_egm.jl:solve_household | Stationary full-grid backward Euler EGM, source floor 1e-5, fine PCHIP | **Different numerical model**: Julia uses 161 active coarse nodes, 641 active fine nodes, no source floor and directly imposed binding-budget policy. Source computes on full 193/769 (including potentially infeasible negative-tail points) and interpolates the constrained coarse policy onto fine grid. Equality not demonstrated. |
| Steady_State_Ayagari/calibrate_R_open_tg.m:16-52 | steady_state.jl:aggregate_prices, stationary_candidate | Factor prices, union labor, GE asset residual | Same economic formula by source inspection; inputs (full vs active EGM/distribution and original eigenvector vs Julia distribution iteration) differ; **numeric equality unverified**. |
| Steady_State_Ayagari/ComputeDistHist_open.m:1-37 | distribution.jl:stationary_distribution | Fine-grid mass splitting and stationary distribution | Both conserve probability and implement Markov transitions; original constructs transition matrix with basefun and eigs, Julia uses iterative distribution; **numerical comparison still required**. |
| Transition/solve_transition_bis.m:7-31 | sudden_stop_equilibrium.jl:stage3_flexible_prices! | Inherited K, labor, production, rental, wages, q, I, profits | Source-level algebra/lag matching largely preserved, but not yet validated at identical numeric candidate K/L. |
| Transition/solve_transition_bis.m:127-153 | sudden_stop_equilibrium.jl:stage3_positive_trade_root | Scalar trade-price root in each period | **Different root-selection protocol** at t=1: source MATLAB fsolve(aggr,0.8), Julia enumerates two positive roots and selects by post-FX occupied household budget feasibility. At observed candidate low 0.7865895867 and high 1.0379166431, these are distinct possible scalar roots. Cross-language branch equivalence unverified. |
| Transition/backsolve_egm.m:22-24 | sudden_stop_equilibrium.jl:stage3_household_backward! | Terminal SS policies for backward induction | **Different dimensions**: source terminal cpol is 193×7 and c_fine 769×7; Julia final steady solution is 161×7 and 641×7 then fills only fine active region. No pointwise full-grid equality yet. |
| Transition/backsolve_egm.m:58-117 | sudden_stop_equilibrium.jl:stage3_household_backward! | Backward EGM given r(t),r(t+1), terms, profits | **Different interpolation**: source EGM on full coarse grid with cutoff Index_b_min, 1e-5 coarse consumption floor, and constrained fine PCHIP; Julia uses active coarse grid, no floor, direct budget at fine constrained points. |
| Transition/backsolve_egm.m:40-55,124-153 | foreign_debt_revaluation.jl:revalue_foreign_debt! | Gross capital, foreign liabilities, exchange valuation and mass splitting | Algebraic source structure preserved (B_shock=0 sets lev=1), with stricter Julia positive-mass grid guards; numerical equality across actual distributions remains unverified. |
| Transition/backsolve_egm.m:168-193 | sudden_stop_equilibrium.jl:stage3_distribution_forward! | Forward mass transport on 769-node fine grid | Same broad Markov and linear mass split, but author's basefun and Julia searchsorted must be matched at boundary/rounding conditions. |
| Transition/solve_transition_bis.m:163-170 | sudden_stop_equilibrium.jl:stage3_market_residual! | 199 asset + 200 union-labor conditions | **Numerical residual difference**: source returns raw market residuals; Julia divides asset by max(1,K_initial) and labor by max(1,L_initial). Zero sets are unchanged by constant nonzero scaling; Jacobian conditioning and numerical values differ. New stage3_unscale_residual! restores source units for comparison. |
| MIT_transition.m:127-154, Main.m:342-351 | foreign_debt_revaluation.jl, scripts/run_sudden_stop.jl | Unleveraged sudden stop before leveraged warm start | Previously wrong order, now corrected; unleveraged solve **still fails**, so leveraged/fixed cannot yet be verified. |
| Transition/solve_transition_fixed.m:5-86 | transition_fixed_exchange.jl | 599 variables incl. 200 terms-of-trade states and Rotemberg Phillips curve | MATLAB includes suspicious last-date price-adjustment signs/indices (source 64-67); Julia follows internally consistent economics rather than necessarily literal MATLAB defects. **Model correction versus numerical reproduction must be kept separate**, and no fixed GE has converged. |

The source and Julia locations above are explicit file/function anchors; exact Julia line numbers may change as implementation evolves. The original MATLAB line numbers are from the GitHub main branch viewed during this audit.

## Confirmed original numerical grid

MATLAB z=linspace(-0.2,1,N) and full asset grid a(z)=bmin+(Amax-bmin) z^3; author explicitly forces z(Index_b_min)=0.
- N=193 coarse: borrowing-limit index 33 (1-based), 32 nodes below borrowing limit, 161 at/above.
- N=769 fine: borrowing-limit index 129 (1-based), 128 nodes below, 641 at/above.
- Julia make_asset_grids already has both complete arrays as coarse_full and fine_full; dropping negative-tail states in the *policy interpolation* is the discrepancy. The integer layout can be checked without solving GE.

**Implication:** It is *not* enough to assign Stage 3's coarse array to coarse_full. Terminal stationary policy dimensions, backward Euler expectation, constrained PCHIP, time-1 FX occupancy, and bound handling must be synchronized. A naive 193-node switch would break the current 161-row terminal policy without establishing source equality.

## Source-reproduction versus economic validity

Source backsolve_egm.m:68-86 floors coarse consumption at 1e-5, while on the fine grid it sets next-period saving to bmin in the binding region. For some negative-tail states these formulas can disagree with the household budget. These states may carry zero or positive mass after FX revaluation. Source behavior and a budget-strict alternative must be reported separately, never mixed silently.

Source scalar trade fsolve from 0.8 does not encode the model selection in Julia's 'prefer budget-feasible root' logic. This is a clear numeric-algorithm difference and may cause non-smooth changes in the GE map as other K/L candidate paths cross its feasibility boundary. A meaningful source-to-source comparison must explicitly report which branch MATLAB returns at each trial.

## New committed controlled reference implementation (no GE equation changes)

- src/source_matlab_audit.jl: source_matlab_grid_layout(...) recovers full grid and Index_b_min; source_matlab_egm_step(...) computes **one** source-style full-grid backward-EGM step given *193×7* terminal policy, with source 1e-5 consumption floor, constrained fine PCHIP and exact borrowing-bound saving. It separately reports fine-grid household budget gaps; it is NOT automatically used by production 399-equation GE.
- src/source_matlab_audit.jl: stage3_unscale_residual!(...) restores author source units from current normalizations, without modifying the equilibrium equations or asserting they are already the same maps.
- matlab/export_source_egm_fixture.m: uses the ORIGINAL repository's rouwenhorst.m to export reproducible full grids, productivity, Markov transitions, one-step source EGM coarse/fine consumption, fine savings and budget gaps for prescribed parameters.
- scripts/compare_matlab_egm_fixture.jl: verifies the exact MATLAB fixture point by point; fails at the first mismatch and prints maximal error and state index.
- test/test_stage3.jl: adds grid index/dimension, source EGM and wrong-terminal-shape tests.

**Important:** These routines are reference/comparison scaffolding, not yet a verified correction to the production GE. The production source household map must not be switched until the reference and MATLAB fixture match.

## Required run sequence and gating criteria

1. Windows Julia syntax/load and unit tests. If they fail, fix the code with tests before proceeding.
2. Run original MATLAB fixture generator from a MATLAB session with the source HANKSOME directory. Do not run Main.m's 399-equation optimization solely to produce a one-date synthetic fixture.
3. Run Julia comparator on the exported MATLAB CSV directory. Report max per-block errors and validate at source-grid index 33 and fine-grid index 129.
4. Replace the production Stage 3 EGM and terminal policies with source-faithful full-grid analogues only **after** the one-date comparison passes and differences to current implementation are measured. Keep source mode separate from budget strict mode.
5. For identical saved K/L paths and identical terminal policy and initial distribution, compare 399 *raw source-unit* GE residuals. If MATLAB's source/Julia diverge, locate the first differing stage, not a 'poor fsolve' symptom.
6. Only after source-level and numerical equivalence, solve zero-leverage 399 unknowns, then 1/16 leverage 399 unknowns, then fixed 599; verify household budgets, full market residuals and figures.

## What is actually verified

The owner verified before this patch: 35 Julia files parse and module load, 2502 unit tests pass, 399-variable unleveraged does **not** converge. The author original source blocks were read and compared by code inspection. A lightweight independent Python count confirms the grid indices implied by linspace, but **neither MATLAB nor Julia are executable in this coding environment**. Therefore the new source audit code, MATLAB fixture, Julia fixture comparator and amended tests are UNEXECUTED and must not be called passed or complete. No new GE computation was attempted.

The central question remains: *Does the Julia 399-dimensional GE residual map match the original MATLAB map on the same complete asset grid and initial conditions?* Answer: not yet; several demonstrable differences remain.

## 2026-11-03: Empirical source verification (supersedes earlier unexecuted note)

The cross-language experiment executes Julia 1.12.6 plus GNU Octave with the ORIGINAL author Exogenous_files/rouwenhorst.m.
GitHub Actions run: https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/actions/runs/41531154396
The synthetic one-step 193/769-grid EGM comparison PASSED. 37 Julia files parsed and loaded, 2420+37+64=2521 tests passed, original shock smoke passed.

| Block | Maximum absolute MATLAB-vs-Julia error |
|---|---:|
| Original coarse grid, 193 | 5.4569682106375694e-12 |
| Original fine grid, 769 | 7.275957614183426e-12 |
| Original 7 productivity levels | 6.719247380715387e-11 |
| Original 7x7 Markov matrix | 6.938893903907228e-18 |
| One-step coarse consumption | 7.257444645247801e-11 |
| One-step fine consumption | 8.558898034749518e-11 |
| One-step fine savings | 1.2005330063402653e-10 |
| One-step fine household budget gaps | 1.2110490388295148e-10 |

The positive asset-increasing synthetic future policy avoids an artificial source PCHIP out-of-support boundary in the earlier constant-consumption probe.

Original-source economic issue: source-consumption flooring and setting constrained fine savings exactly to bmin can produce large negative household budget gaps on below-bmin current asset nodes (around -197 in the synthetic probe). Those points may be unoccupied; positive mass after FX revaluation must be checked before claiming a valid economic equilibrium. Original-source numerical matching and a budget-valid model must be reported separately.

Stationary reference: Julia and Octave each iterated the original full-grid solve_R-style policy to convergence after 933 iterations for synthetic parameters. The coarse stationary policy matched to 7.257444645247801e-11. The original Octave stationary fine-grid c_fine contained FOUR NA states. The first Julia reference erroneously extrapolated at these states rather than preserving original interp1 undefined output. A follow-up patch returns source NaN outside the endogenous knot range, and comparator reports NaN mask equality; that second check is pending.

DO NOT infer full calibrated steady-state or 399-variable GE matching. These are controlled source one-step and synthetic stationary experiments, not the author's full calibrated source simulation.

Remaining gates: original stationary full 193/769 policy under author prices, production Stage3 terminal policy matching, FX root branch from source fsolve(0.8), distribution and raw 399 market residual equality. Do not run another 399 root solver until these gates are established.
PR #9 remains Draft.


## Exact run evidence, 2026-11-03 — full-grid error is too small to explain initial 399 residual

Read CI job log of [run 41539508608](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/actions/runs/41539508608). Both Julia 1.12.6 and the author's MATLAB routines under GNU Octave ran successfully. The full source-style calibrated policy (193 coarse, 769 fine) **matched the original solve_R.m** at calibrated Julia prices: maximum coarse error 1.1368683772161603e-12 and maximum finite fine error 1.1937117960769683e-12. Comparing the source policy against the existing **161-node production** stationary policy, the active coarse maximum difference was only 2.984211278089788e-8; active fine maximum 2.4454461765632818e-5.

On the same **initial candidate** K/L path, the entire full-grid source-style backward transition and independently propagated full-grid 769-node distribution produced:
- Maximum **scaled asset residual-map difference**, full source-style minus production: **5.0586683097967915e-9** (period 199).
- Maximum **scaled labor residual-map difference**: **4.7197567187851064e-8** (period 1).
- Period-1 aggregate assets: source 21.665907683451184 vs production 21.665907578130664.
- Period-1 domestic consumption: source 1.3952003883065274 vs production 1.3952004514988403.
- Maximum all-node fine policy difference around 0.1948–0.2412, but the difference over **occupied household nodes** was at most about 2.58e-5 and their mass-weighted difference about 6.3e-8 at period 1. Zero-mass tail policy differences must not be confused with their effect on equilibrium residuals.
- Occupied source-invalid mass **0**, maximum occupied budget gap **<=2.5e-16** in the tested candidate path.

**Decision:** The 161-versus-193 coarse grid is a REAL source algorithm difference, but at the tested initial path its effect on the 399 residual map is many orders of magnitude smaller than the observed minimum unsolved residual (~0.005). Do **not** replace the production grid or interpret it as the proven root cause of nonconvergence on this evidence. A source-map comparison at the **last accepted but unconverged checkpoint** is the next discriminating observation; the original-versus-production difference may grow along the GE iterations.

**Original author trade-price branch:** At the same opening candidate, the unmodified original scalar `fsolve(aggr,0.8)` returned **0.7865895853334468** for both no-leverage and leverage=1/16. The analytic roots are ~0.7865895866917804 and 1.037916643090368. Current budget-strict Julia selects the lower root with no leverage but **the upper root under leverage=1/16**, because the lower root puts positive-mass households into an infeasible budget state (min resources -0.1826131304704792 in owner log). This is a **verified source branch mismatch**, distinct from the grid difference. The original source's consumption floor cannot automatically be called a valid economic solution.

**Original basefun FX map numerical equality:** For the real shock-date 769x7 distribution, maximum distribution differences between the original MATLAB basefun and Julia were 3.56e-17 (no leverage low), 9.97e-18 (no leverage high), 7.98e-17 (1/16 low) and 1.21e-17 (1/16 high). The FX mass-splitting implementation is therefore validated for these candidate roots; generic revaluation debugging is not the immediate bottleneck.

This successful experiment evaluated source-faithful **Julia reference** residual maps, not yet the original `solve_transition_bis.m` 399-vector end-to-end. Direct source output can still diverge because of 0*NaN propagation and fsolve branch details. The next commit adds a direct one-pass source `solve_transition_bis.m` evaluation using identical Julia-exported 399-vector, B path, terminal 193-node policy and household distribution. It must FAIL loudly for nonfinite or mismatching source residuals, rather than implying full replication.

**Next discriminating command on owner's Windows system (without GE root finding):**

```powershell
julia --project=. scripts/export_stage3_fx_fixture.jl results/checkpoint_399_audit --checkpoint
```

This reads `results/stage3_unleveraged_UNCONVERGED.checkpoint`, recomputes the current Julia production map and original-source-like full-grid EGM/distribution at that **same 399-vector**, and exports per-period residual discrepancies. This checkpoint is an accepted iterate, **not** an equilibrium. The source-vs-production 399-map difference there will determine whether the full-grid discretization materially explains the actual nonlinear solver failure.
