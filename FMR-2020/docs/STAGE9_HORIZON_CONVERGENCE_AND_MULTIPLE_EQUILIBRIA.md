# HANKSOME Stage 9 — 200-to-300 horizon convergence and multiplicity certification

**Status:** The 200-period strict-portfolio appreciation and depreciation **market** roots are independently verified. A 300-period **899-equation version of the IDENTICAL initial household economy and foreign-credit shock** has been constructed and both prolonged Stage 8 roots have been successfully evaluated with the current Euler-polished household and strict portfolio GE residual. **No T=300 nonlinear market root has yet been solved.** Infinite-horizon multiplicity is not certified.

Target: Sergio de Ferra, Kurt Mitman and Federica Romei (2020), *Household Heterogeneity and the Transmission of Foreign Shocks*, *Journal of International Economics* 124, 103303, https://doi.org/10.1016/j.jinteco.2020.103303.

Previous evidence: `docs/STAGE8_TWO_MARKET_ROOTS_AND_TERMINAL_CERTIFICATION.md`; `results/stage8_depreciation/two_market_root_comparison.csv`.

## 1. Why horizon extension is needed

The **two opposite-sign market roots at T=200 belong to the same initial inherited state and shock** when re-evaluated with the existing strict-portfolio Euler-polished household GE. They are numerically different, but the supplied date-200 stationary continuation is not attained perfectly.

| T=200 terminal / market diagnostic | Appreciation | Depreciation |
|---|---:|---:|
| Impact terms x1 | 1.156101137755 | 0.868649787341 |
| 599 max scaled residual | 1.2302851e-7 | 3.2355748e-7 |
| Nominal depreciation (log points) | −11.82727 | +16.75958 |
| Terminal household budget max on all positive mass | 0.011342185 | 0.076587683 |
| Terminal budget max on states with mass >1e−12 | 0.002880380 | 0.020738305 |
| Probability mass with terminal budget error >1e−4 | 70.61319% | 71.12532% |
| Probability mass with terminal budget error >1e−3 | 0.30702% | **26.78114%** |
| Terminal mass-weighted absolute household budget error | 0.000180931 | 0.000824717 |
| Max stationary copied-policy budget error | 4.55e−13 | 4.55e−13 |

The worst all-positive-mass terminal node has probability mass 5e−324, but **material mass states also violate the terminal continuation**, especially along the depreciation branch. Thus excluding an almost-zero mass node alone cannot fix economic certification. Neither the distribution nor its moments are truncated.

The Stage 8 identity for each terminal household is exactly

```text
R_T = R_terminal_stationary_policy
    + (r_T-r_ss)*a
    + (w_T*L_T-w_ss*L_ss)*s
    + (profits_T-profits_ss)
    - (m(x_T)-m(x_ss))*c_terminal_stationary.
```

The stationary-policy budget error is effectively zero. At T=200, the actual-vs-stationary interest gap (especially at large a) is the main source of the depreciation branch terminal discrepancy. We must test **whether this and the material-mass gap decrease when T is raised**.

## 2. Exact same economy, only its numerical truncation changes

We reuse the already implemented and tested `FMR2020.stage3_extend_horizon` and `FMR2020.stage3_prolong_aggregate_guess`; no new GE, EGM or distribution algorithm is introduced.

The inherited state is reconstructed **once** from the original saved T=200 Euler-polished Julia credit-expansion root using `stage5_existing_boom` and `stage3_reversal_context(steady, boom; shock_start=41, leverage_hat=1/16)`. This is the same origin as both Stage 8 market roots, not the different original MATLAB source state.

The 300-period context:
- Is identical at ALL T=200 dates in aggregate exogenous foreign bond supply `B_t` and aggregate productivity `z_t`.
- Holds `B_{201},...,B_{300}=B_{terminal}` exactly. The original source shock reaches its calibrated terminal credit supply at T=200; **it is not stretched/redefined to reach steady state only at T=300**.
- Retains exactly the original household initial assets and productivity distribution, borrowing limit, full fine and coarse grids, transition matrix, gross portfolio scaling, initial capital, inherited B and q, and pre-shock terms.
- Retains exactly the same calibrated terminal steady state, including its household terminal policy.
- Has `899=299+299+299+2` equations: asset/labor/actual realized home-goods markets for dates 1–299, plus terminal labor and terms anchors.

The original 200-date K/L/terms paths are copied *unchanged* into their 300-date initial guesses; new K is geometrically connected to the same terminal stationary K, with new L and terms initialized to their stationary levels. This **prolongation is only a numerical starting point**, not a purported market root.

All of these invariants are checked in `test/test_stage9_horizon.jl` and again before starting the real 899 nonlinear solve.

## 3. Actual first 899-vector residual evaluations (no solve)

The Julia 1.12.6 Stage 9 Actions audit evaluated both prolonged root guesses under the **current full 899 residual**, including all 300 Euler-polished household policies, foreign-debt revaluation and forward distributions. In both cases the full residual map was defined.

| T=300 extended initial candidate | Appreciation-side | Depreciation-side |
|---|---:|---:|
| Impact terms of trade x1 | 1.156101137755 | 0.868649787341 |
| Maximum scaled 899 residual | **0.001155321664** | **0.000378689454** |
| Maximum scaled asset residual | 0.000550764312 | 0.000171737888 |
| Maximum scaled labor residual | 0.001155321664 | 0.000378689454 |
| Maximum scaled goods residual | 0.000607963085 | 0.000217895275 |
| Actual first-period scaled goods residual | −3.94050e−8 | +1.37185e−8 |
| Interior maximum occupied household budget | 4.44e−16 | 4.44e−16 |
| Probability mass conservation max error | 1.33e−14 | 1.32e−14 |
| Terminal budget maximum over ALL positive mass | 4.77775e−5 | 1.79239e−5 |
| Terminal budget maximum on cell mass >1e−12 | 1.26876e−5 | 4.86028e−6 |
| Mass-weighted absolute terminal budget | 4.54565e−7 | 1.74081e−7 |
| Terminal budget violation mass >1e−4 | **0%** | **0%** |
| Terminal budget violation mass >1e−3 | **0%** | **0%** |
| Terminal interest-rate deviation | +1.76501e−8 | −6.76128e−9 |
| Terminal capital gap | +2.77732e−7 | −1.06392e−7 |
| Status | **Defined, NOT a T300 root** | **Defined, NOT a T300 root** |

Direct source: successful [Julia CI run 41639822611](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/actions/runs/41639822611). The actual 18-field per-branch values are checked in as [verified_initial_candidate_summary.csv](../results/stage9_horizon/verified_initial_candidate_summary.csv), transcribed from the CI console output with no new solve. Complete metric/value files are archived as artifacts of the same CI run. **No 899-dimensional nonlinear root solving was carried out in CI.**

The extended numerical initializations have already reduced the **weighted absolute** terminal budget errors by approximately 398 times (appreciation) and 4,738 times (depreciation), *before solving any 899-equation GE*. This is evidence consistent with a truncation issue rather than a defective stationary policy. It is **not yet a horizon-convergence finding about T300 equilibria**; the full 899-market residual must be reduced below tolerance and the resulting path independently rechecked.

The corresponding source code is `scripts/stage9_horizon_extension.jl`. The original 599-vector root CSV files are never rewritten.

## 4. Validated implementation and numerical work still required

The Stage 9 driver provides FOUR explicit modes:

```powershell
Set-Location C:\github\22p8z0omwz9xu
git switch research/fmr-stage9-horizon-multiplicity-certification
git pull --ff-only
Set-Location FMR-2020
julia --project=. scripts/check_syntax.jl
julia --project=. test/test_stage9_horizon.jl
julia --project=. scripts/stage9_horizon_extension.jl --audit
julia --project=. scripts/stage9_horizon_extension.jl --solve-depreciation
julia --project=. scripts/stage9_horizon_extension.jl --solve-appreciation
julia --project=. scripts/stage9_horizon_extension.jl --compare
```

The **default and `--audit`** mode evaluates both 899-vector starting candidates, saving (and CI archiving) terminal mass-aware and current market residuals. The owner-PC **`--solve-depreciation` and `--solve-appreciation`** modes use **the existing `solve_stage3_flexible_full`** with `strict_portfolio=true`, `polish_households=true`, tolerance 1e−6, `maxiter=70`, `max_evaluations=5000`. No artificial five-step cutoff or altered initial shock is introduced.

Every accepted solve step is stored as a separate 899-element file plus a restartable checkpoint; the log records impact terms, currency response sign, full residual and worst equation. If an interrupted PC attempt is resumed, the saved residual is first recomputed with the CURRENT 899 map to avoid stale-model restarts. Convergence to the *opposite* currency-sign root is not silently relabeled.

A successfully solved 899 market root is exported under `*_T300_MARKET_ROOT_UNCERTIFIED.csv`, with independent household-budget, Euler/KKT, probability mass, asset/labor/goods market and terminal diagnostics in `*_T300_economic_certification.csv`. The filename explicitly remains **UNCERTIFIED** until terminal/horizon properties have also been verified.

The **`--compare`** mode only runs once both 899 roots exist. It re-evaluates old and new roots under their respective unchanged economies, reports the two impact terms, 200-vs-300 terminal budgets by mass, and common first-40/100/199 calendar-date K/L/terms, home consumption, output and investment differences.

New heavy 899-variable GE runs are intentionally **not** performed on GitHub Actions; they must be performed on the owner's PC.

## 5. Acceptance decisions and unresolved economic questions

The Stage 9 results must establish more than new market roots. Compare interest rate, capital and wage convergence at the end of 300 dates, occupied household consumption/borrowing/Euler/KKT, market clearing and terminal budget violations over all positive mass and on material probability mass support.

- **Outcome A:** Both opposite-price market roots survive T300, share precisely the same original economy, and terminal/horizon diagnostics show sufficient convergence for strong (but not theorem-level) evidence of infinite-horizon multiplicity.
- **Outcome B:** Both T300 market roots survive but terminal mismatch remains economically material.
- **Outcome C:** Only one T300 root is numerically certified; failure to solve the other is not a proof of nonexistence.
- **Outcome D:** The impact and path allocations converge toward the same branch, indicating truncation sensitivity.
- **Outcome E:** Neither branch can be certified or numerical feasibility prevents the required comparison.

**Current status before the owner-PC 899 solves: Outcome E (full T300 certification unresolved).** The original T200 two-root result is established, and both 899 starting vectors evaluate successfully. Do not claim A/B/C/D without actual T300 numerical solutions.

The original paper studies the welfare and heterogeneous foreign-currency balance-sheet implications of depreciations and capital-flow reversals, but its descriptive transmission results by themselves do not establish uniqueness or prescribe which of two newly discovered strict-portfolio roots the private sector selects. A numerical starting-value basin is not an expectations-selection mechanism. A definitive statement about the original paper's uniqueness treatment requires direct verification of its equilibrium appendix; no new uniqueness or determinacy theorem is claimed here.

## 6. Code and validation provenance

- Branch: `research/fmr-stage9-horizon-multiplicity-certification`.
- T200 source roots: `results/economic_reproduction/economic_flexible_leveraged_euler_polished_STRICT_PORTFOLIO_MARKET_ROOT_UNCERTIFIED.csv` and `results/stage8_depreciation/source_initialized_STRICT_MARKET_ROOT_UNCERTIFIED.csv`.
- Implementation: `scripts/stage9_horizon_extension.jl`, reusing existing Julia 3T−1 solver and Stage 8 terminal household audit.
- Regression tests: `test/test_stage9_horizon.jl`, verifying horizon dimensions, unchanged inherited distribution and full B path, exact old-path prefix, two root signs, and explicit dimension failures.
- Workflow: `.github/workflows/fmr-stage9-horizon-audit.yml`; Julia 1.12.6, `Meta.parseall`, module load, exact state/horizon tests, and actual 899 residual evaluation ONLY.
- No original MATLAB .mat files or Stage 5–8 equilibrium CSV/checkpoints are overwritten.

At the completion of actual T300 solves, update this report with both new roots and all independent terminal/horizon numbers. If only an unverified 899 root is obtained, report it as a provisional market root, not a complete infinite-horizon equilibrium.
