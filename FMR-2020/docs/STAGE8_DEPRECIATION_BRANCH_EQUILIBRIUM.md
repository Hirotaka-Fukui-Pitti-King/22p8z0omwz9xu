# HANKSOME Stage 8 — Depreciation-branch GE: original-path strict-system audit and controlled numerical decision

**Research verdict (as of 2026-11-03): Outcome E. The original MATLAB depreciation candidate is an economically feasible starting point for the corrected strict-portfolio Julia problem but is not a full-market root. A separate valid depreciation equilibrium is not established or excluded. No new 599-variable nonlinear solution was run in the CI audit.**

Paper: de Ferra, Mitman, and Romei (2020), *Journal of International Economics* 124, 103303, https://doi.org/10.1016/j.jinteco.2020.103303. Original source: https://github.com/kurtmitman/HANKSOME. Working scientific context: Stage 7 `STAGE7_CONSTRAINED_POLICY_AND_EQUILIBRIUM_RECONCILIATION.md`.

## 0. Owner-PC Stage 8 full 599 solve and checkpoint recovery (2026-11-03)

The user has now **actually executed** `julia --project=. scripts/stage8_depreciation_continuation.jl --solve` on Windows 11, Julia 1.12.6, using the original depreciation path as the starting vector, with exactly the documented strict-portfolio and Euler-polished full 599-equation map.

The solver accepted **five** steps, with **no full Jacobian rebuild failure and no exhausted function-evaluation budget**. Its `maxiter=5` stopped the run by design:

| Accepted full-system iteration | Maximum scaled 599 residual | Damping |
|---:|---:|---:|
| 0 | 0.107256061729 | — |
| 1 | 0.077513986653 | 0.25 |
| 2 | 0.025723462623 | 1.0 |
| 3 | 0.005446066350 | 1.0 |
| 4 | 0.001638584342 | 1.0 |
| 5 | **0.001520252578** | 1.0 |

The last solver call reports `evaluations=609`, while its callback has already saved the fifth accepted **full 599-vector iterate** as `results/stage8_depreciation/source_initialized_UNCONVERGED.checkpoint`. This is **not** a failed feasibility check, a nonexistence result, or a converged equilibrium. It is an interrupted but improving full-system Newton/Broyden sequence.

**New recovery implementation:** The Stage 8 script accepts explicit `--resume`. It uses `stage3_read_checkpoint` (read-only), independently reconstructs the same strict inherited Julia boom and the corrected household/market operator, recomputes all 599 residuals, and **refuses stale-model checkpoints** unless the recorded and actual maximum residual agree to `atol=1e-8, rtol=1e-5`. It skips the repeated source 3×3 probe and begins a bounded **12-step, at most 2200 residual-evaluation** follow-on call to the existing solver from the actual saved fifth iterate, not the original MATLAB guess. **It reinitializes the 599-column Jacobian**, because the existing checkpoint contains x only and does not save the previous Broyden matrix. It does not claim continuous Broyden state.

Each new accepted iterate records the current impact trade price, nominal depreciation log points, maximum *full* residual, residual norm and most binding equation in `results/stage8_depreciation/accepted_full599_trace.csv`. Before resuming, the originally accepted iteration-5 full vector is copied unchanged to `accepted_iteration_0005.checkpoint`; every subsequent accepted vector is also stored separately under `accepted_iteration_NNNN.checkpoint` before the rolling resume checkpoint is replaced. This preserves numerical branch history, even if the later solver approaches appreciation (original MAT and Stage 5 root remain untouched). If the iterate crosses into appreciation, the trace labels that transition; convergence to the already-known appreciation solution would **not** establish a distinct depreciation equilibrium.

On the owner PC, from the repository root:

```powershell
git switch research/fmr-stage8-depreciation-continuation
git pull --ff-only
Set-Location FMR-2020
Get-Content .\results\stage8_depreciation\source_initialized_UNCONVERGED.checkpoint -TotalCount 4
julia --project=. scripts\stage8_depreciation_continuation.jl --resume
```

The four header lines must begin `# STAGE3_UNCONVERGED_CHECKPOINT v1`, `# n=599`, `# iter=5` and `# max_scaled_residual=0.001520252578...`. The explicit `--solve` source restart now REFUSES to overwrite an existing checkpoint, in order to prevent accidental loss of the accepted five-step progress.

**Current actual classification remains Outcome E**, pending the first resume. Do not claim that `--resume` has been run on the user's PC or that a low-price root has converged. Only the read/validation and unit-test path can be checked in lightweight CI.

## 1. Precisely which economy was evaluated

The **target economy** is the current Julia *Euler-polished, strict gross-portfolio, flexible-price, inherited surprise-shock* general equilibrium:

- The inherited pre-shock full asset/productivity distribution is **reconstructed from the actual archived 200-period Julia Euler-polished credit-expansion aggregate GE root** (not the original MATLAB saved distribution). The archived boom root is re-evaluated once, not nonlinear-solved again; its full endogenous pre-shock distribution is regenerated from household EGM and forward transport.
- At unexpected credit reversal date 41, inherited `K_1`, pre-shock shadow price `q_0`, foreign bonds `B_0`, the distribution, and pre-shock terms `x_pre` all come from that same saved Julia boom, using `stage3_reversal_context`.
- Revaluation uses `foreign_bonds*x_pre/x_1`. **`strict_portfolio=true`** rescales household gross capital to match the inherited `q_0 K_1`, while household bonds exactly account for the remaining net wealth and inherited external bonds.
- The 200-period exogenous `B_t` is the exact already-implemented quartic credit reversal back to the calibrated terminal foreign-asset supply. No residual-based or hidden parameter interpolation is used.
- The household map is the **existing** `stage3_household_backward!` and `stage3_euler_polish_households!`, with no consumption floors for positive-mass states. Actual household consumption from the **same** distribution is used inside every interior goods-market residual.
- The candidate is `X=(K_2,\ldots,K_{200};L_1,\ldots,L_{200};x_1,\ldots,x_{200})`, with **199 asset, 199 labor, 199 actual household home-goods residuals, and two terminal anchors**. Date-200 asset/goods equations are not falsely imposed against the exogenous stationary consumption boundary.
- First-period firm prices and investment costs are recomputed from this strict-economy `K_1` and the candidate `K_2,L_1,x_1`. No original MATLAB firm price/output is transplanted.

**Important difference of initial conditions:** the original MATLAB source `K_1=29.75127255693537`, `q_0=1.0171493435307144`, and `x_pre=1.0258494514359664` are *not* the same as the Julia archived boom `K_1=29.593930147144125`, `q_0=1.0172430715493037`, and `x_pre=1.0271423614420732`. Both sources use `B_0=-6.281005741472044` in this comparison. The MATLAB `K_2:L_{200}:x_{200}` vectors are therefore **numerical initial guesses only**, and not a distinct equilibrium of the Julia initial-state economy.

The independently solved Stage 5 strict-portfolio appreciation root, reported in previous PC logs as `x_1=1.1561011378`, was **not recomputed**. Its strict-root CSV is not currently tracked in this repository, so a same-state comparison with its entire 599-vector cannot be independently certified on CI until that archived file is supplied. Legacy portfolio market-root files are *not* substituted for it.

## 2. Original MATLAB depreciation candidate under the strict 599 map

Using `scripts/stage8_depreciation_continuation.jl` in default audit mode, Julia 1.12.6 GitHub Actions reconstructed the saved Julia pre-shock GE, imported the *actual original* `leveraged_candidate.csv` and `leveraged_source_terms.csv`, and computed one **strict-portfolio, Euler-polished** 599-vector. The source estimate `x_1=0.90985851141107 < x_pre=1.02714236144207` is a genuine depreciation numerical guess.

| Source-initialized strict-portfolio candidate | Verified CI result |
|---|---:|
| Impact source terms `x_1` | 0.909858511411 |
| Original candidate `K_2` | 29.672235297975 |
| Original candidate `L_1` | 1.065590530478 |
| Current inherited foreign B | −6.281005741472 |
| Terminal foreign B | −4.343249891001 |
| Maximum scaled asset-market residual | **0.0317140906994** (date 1) |
| Maximum scaled labor-market residual | **0.1072560617286** (date 1) |
| Maximum scaled actual goods-market residual | **0.0913646252006** (date 1) |
| Maximum scaled 599-vector residual | **0.1072560617286** |
| First-period goods residual in model units | **−0.279808550185** |
| First-period actual household home-good demand | 1.451794780937 |
| Scaled terminal labor anchor | −4.92221731363e−6 |
| Scaled terminal trade-price anchor | −0.00123504454691 |
| Maximum distribution-mass error | **1.24344978758e−14** |
| Maximum positive-mass interior household budget residual | **4.4408920985e−16** |
| Minimum positive-mass interior consumption | **+0.157901735422** |
| Positive-mass interior borrowing/consumption violations | **0%** |
| Date-200 capital gap against terminal SS | **−0.006571897459** |
| Minimum original-state household resources *before* wealth-grid interpolation | **+0.2736715795904** |
| Probability mass with nonpositive original-state household resources | **0%** |
| Probability mass revalued outside asset-grid support | **0%** |
| Strict portfolio rescaling relative to raw source-style gross capital | 1.0010605649503 |
| Raw gross-capital accounting gap at this initial state | −0.0318936557216 |

**Key conclusion:** the source-depreciation candidate is not rejected by the shock-date *necessary original-state household affordability test*, and the entire polished 599 residual can be evaluated. But asset, labor and goods markets each fail by economically meaningful amounts. It is a defined **initial guess**, not a solved GE.

The raw portfolio gap shown above is deliberately measured **before** applying the correct scaling; it is not an uncorrected gross-capital failure in the strict system.

## 3. First-date full-household-map Jacobian and a single damped local step

The source candidate's actual 599 residual was differentiated at three aggregate coordinates **only**, with the other 596 coordinates held fixed. Each finite-difference probe reruns the entire endogenous household policy, shock revaluation, distribution transport, and goods/labor/asset residual map with `polish_households=true`, `strict_portfolio=true`.

The computed impact goods-market derivatives are

```text
d R_H(1) / d K_2 = −0.5884040333310
d R_H(1) / d L_1 = +0.6209174433099
d R_H(1) / d x_1 = +1.316382920677
```

The full local 3×3 Jacobian has rows `(asset_1, labor_1, goods_1)` and columns `(K_2,L_1,x_1)`:

```text
             K2             L1             x1
asset_1   -0.184012938    +0.054686738    +0.239457215
labor_1   -0.999096612    -1.732230409    -0.543052226
goods_1   -0.588404033    +0.620917443    +1.316382921
```

Its condition number is approximately `125.9887407`. The local Newton direction from the original source guess is

```text
delta K2 = -0.3351972381790
delta L1 = +0.1838056015856
delta x1 = -0.1671207385423.
```

The undamped step attempts further depreciation and fails an occupied-household feasibility check: the worst ORIGINAL pre-interpolation state has disposable resources **−0.1269082063905**, with violating probability mass **0.0102476895031 (1.02477%)**, even though its revalued wealth remains within the saved asset grid. It is **not** an admissible economic solution. With damping `alpha=0.5`, the full household/residual map is finite, and the accepted **local three-equation** numerical trial is

```text
K2 = 29.504636678886
L1 = 1.157493331271
x1 = 0.826298142140  (still depreciation).
```

The damped trial's minimum original-state disposable resources are **+0.0667983222232**, with zero violating original-state mass. The *norm of the first three residuals* changes from `0.144420016149` to `0.081740984141`. The trial's first-date asset/labor/goods residuals are `(-0.0147676877201,-0.0273816513863,-0.0755893448520)`.

**However, the maximum scaled FULL 599 residual rises from `0.107256061729` to `0.131528476142`.** This trial is therefore **not** a globalization step for the actual 599-variable general equilibrium, let alone a root. Improving the first three equations while holding all future aggregates fixed can worsen the remaining dynamic markets. This observation demonstrates why solving a scalar price root or a three-equation first-period problem cannot certify the full depreciation branch.

These are local numerical derivatives on a fixed inherited state, **not** general-equilibrium comparative statics or a statement about equilibrium selection.

## 4. Household budget and terminal tests

The script separately audits interior occupied household budgets, borrowing minima, and distribution mass using the retained `TransitionWorkspace` from the strict residual computation. These identities are independent of the solver convergence flag; in particular, they do not convert an off-equilibrium source guess into a fictitious converged `TransitionSolveResult`.

The initial-source terminal labor and terms anchors are near their targets, but this does not certify terminal stationary household budgets, full terminal capital consistency or infinite-horizon truncation. Euler-polished policies are recomputed inside every candidate residual; the script does **not** claim an independent all-period Euler/KKT certification or the existence of a full market root at the source guess.

## 5. Computational decision and owner-PC command

The source-initialized candidate passes necessary household affordability, so a **single bounded 599-variable solve** from this actual low-price guess is economically justified. That solve is deliberately **not run in GitHub CI**: the project explicitly reserves expensive 599-dimensional finite-difference Jacobian and Newton computations for the owner's Windows PC. This ChatGPT environment has no access to that user's local Julia runtime.

The first five-step full-system attempt has been executed on the owner PC and ended at the preset iteration limit, as documented in §0. The source-initialized first-run command remains reproducible, but **must not be repeated when its checkpoint exists**:

```powershell
Set-Location C:\github\22p8z0omwz9xu
git fetch origin
git switch research/fmr-stage8-depreciation-continuation
git pull --ff-only
Set-Location FMR-2020
julia --project=. scripts/check_syntax.jl
julia --project=. test/test_stage8_depreciation.jl
julia --project=. scripts/stage8_depreciation_continuation.jl
julia --project=. scripts/stage8_depreciation_continuation.jl --solve
```

The **default** command is evaluation-only. The original `--solve` command (already executed once on the owner PC) calls the **existing** `solve_stage3_flexible_full`, with `polish_households=true`, `strict_portfolio=true`, tolerance `1e-6`, maximum **5** Newton/Broyden iterations and **1000** full residual evaluations. Accepted-but-unconverged iterates are checkpointed separately under `results/stage8_depreciation/source_initialized_UNCONVERGED.checkpoint`; no prior Stage 5 root or original MATLAB file is replaced.

If the bounded solver converges, `stage3_diagnostics` independently checks household budgets, market levels and probability mass before any new market-root CSV is written. Its filename explicitly says `MARKET_ROOT_UNCERTIFIED`. A converged finite-horizon market map is still **not** evidence of a certified infinite-horizon equilibrium; positive-mass Euler/KKT, terminal policy/budget and horizon sensitivity must then be audited.

If the bounded solve fails, record the checkpoint, last accepted scaled residual, worst market, and feasibility diagnostics rather than rerunning an unrestricted multistart. Do not automatically alter the foreign-debt shock or the household borrowing constraint. Any economically consistent foreign-credit homotopy requires matching the initial inherited B and terminal steady state at each intermediate path; no ad hoc interpolated residual system or purported continuation result is claimed by the current code.

## 6. Result classification

**Outcome E — No newly certified full depreciation-side equilibrium.** A single original-source numerical guess in the strict Julia economy was evaluated and its first-date Jacobian analyzed. The source guess is financially feasible at positive-mass original shock-date states, but the full 599-equation residual is large. A half-step can reduce the three impact market errors while preserving depreciation, but **increases** the full dynamic residual. An undamped more-depreciated trial is infeasible on a positive-mass state. These facts neither prove that depreciation equilibrium does not exist, nor demonstrate a second equilibrium.

The already known Stage 5 strict appreciation path and the original MATLAB source depreciation path cannot be called two equilibria of the same model without independently solving the 599-system for the *same* corrected initial state and checking every market, household and terminal condition.

**Remaining crux:** can the existing full-Jacobian 599 nonlinear solver jointly adjust capital, labor and the entire trade-price path from the now verified feasible low-price initial guess to a second **financially feasible** market root? The saved evidence narrows the question to full dynamic GE, not the Stage 7 interpolation error or the source nested scalar minimizer.

## 7. Recorded artifacts and verification

- `scripts/stage8_depreciation_continuation.jl`: original-source numerical initialization, strict-portfolio reference economy reconstruction, one full 599 evaluation, local sensitivity and bounded owner-PC-only solve.
- `test/test_stage8_depreciation.jl`: source 399/200/599 dimensions, archived actual MATLAB first price, serialization and non-certification.
- `.github/workflows/fmr-stage8-depreciation-audit.yml`: Julia 1.12.6, package instantiation, full source syntax parse, test, **one** strict source-initialized full residual evaluation and optional at-most-small first-date Jacobian probes; **no full GE nonlinear solve**.
- The tracked `results/stage8_depreciation/initial_candidate_residuals.csv` contains 90 observed metric/value pairs, copied exactly from the successful [Julia 1.12.6 audit run 41625705198](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/actions/runs/41625705198). The CI artifact `fmr-stage8-strict-source-candidate-audit` also includes the direct script output. No GE root or numerical estimate was fabricated.
- Existing original author `.mat`, saved Stage 5 roots and prior Stage 6/7 audits remain intact.

### Important verification qualifications

Julia syntax parsing and a successful audit run on Julia 1.12.6 verify the **default evaluation mode**. The opt-in **`--solve` numerical branch WAS executed on the owner PC for five accepted steps and stopped at maxiter=5 with residual 0.001520252578, without convergence**. The newly added **`--resume`** numerical branch has not been executed on the owner's PC at the time of this update. Light CI checks do not constitute a resumed full solve. The preexisting stricter Euler refinement and full horizon convergence remain questions for any future full root.
