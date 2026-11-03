# Stage 6 — Independent equilibrium audit of the ORIGINAL MATLAB depreciation path

## Decision and reference

This document audits the original authors' **saved** leveraged flexible-FX sudden-stop experiment (`FMR-2020/results/author_reference/exp_flex_lev16.mat`) from de Ferra, Mitman and Romei (2020), *Journal of International Economics* 124, 103303. The primary code-level reference is the **original** `HANKSOME/Transition/solve_transition_bis.m` and `HANKSOME/Transition/backsolve_egm.m`.

The original MATLAB depreciation path **does not satisfy the independent first-period goods-market or full first-period household budget identities to numerical equilibrium accuracy**. The original household balance sheets are, however, **financially feasible** at the authors' own saved prices, wages, interest rates and profits. This separates the actual MATLAB path from the different *conditional* low-price root previously examined with the Julia appreciation-branch aggregate allocation.

These numerical findings reject the stored path as a *certified full economic equilibrium at its saved numerical precision*. They do **not** prove that a corrected depreciation-side full GE is impossible, nor that every original economic mechanism is wrong. Identification of the underlying source of the residuals, whether a genuinely distinct all-market low-price GE exists, and infinite-horizon terminal validity remain separate questions.

## Original-source calculations, no MATLAB/Julia GE solve

GitHub Actions [Stage 6 audit](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/actions/workflows/fmr-stage6-author-audit.yml) runs `scripts/stage6_author_depreciation_audit.py` with Python, NumPy and SciPy. The audit opens the tracked original `.mat` file directly using `scipy.io.loadmat`, rather than reconstructing a synthetic economy or calling MATLAB's `fsolve`. Outputs are attached to each CI run as `fmr-stage6-original-matlab-depreciation-audit`, containing `first_period_metrics.csv`, `first_period_occupied_budget_errors.csv`, `first_30_periods.csv` and `stage6_audit.md`. The routine does **not** change the author's saved data, any portfolio mapping, the MATLAB code, or the archived Julia market roots.

Audited experiment: flexible FX, capital-intercept 1/16, shock immediately after period 40 of the previously anticipated credit expansion, saved 200-period reversal.

### Observed results — original MATLAB saved allocation

| Diagnostic | Actual value |
|---|---:|
| Before-shock terms of trade `x_pre` | 1.025849451436 |
| Original shock-date terms `x_1` | 0.909858511411 |
| Nominal depreciation, `100 log(x_pre/x_1)` | **+11.99871760 log points** |
| Maximum original 399-equation outer residual (unscaled) | 1.3536933352e-8 |
| Source scalar `c_t_guess(x_1)-net_cons(x_1)` | **+0.00211041927838** |
| Exact source-scalar minimizing positive price | 0.909875606122 |
| Exact global minimum of source-scalar residual | **+0.00211041785147** |
| Source-scalar positive zero at the saved macro path | **Does not exist** |
| Actual household consumption minus source aggregate-budget consumption | **+0.00222480837002** |
| Actual first-period household consumption `C_H` | 1.33608285423 |
| Actual home-goods residual `net_cons(x_1)-C_H` | **−0.00433522764840** |
| First-period asset-market residual, level | +1.29255752768e-8 |
| First-period labor residual (source outer equation), level | −5.39561950674e-10 |
| Gross capital holdings minus required `q_0 K_0` | **−0.0318889358579** |
| Gross foreign bonds minus inherited `B_0` | **+0.0318889491319** |
| Aggregate initial net wealth residual | +1.32740076708e-8 |
| Lowest feasible original-state resources, before interpolation | **+0.243344083414** |
| Probability mass violating necessary affordability before interpolation | **0** |
| Off-grid probability mass after wealth revaluation | 0 |
| Lowest post-interpolation grid-state resources | +0.232172443672 |
| Post-interpolation probability mass violating affordability | 0 |
| First-period maximum individual budget residual (absolute) | **0.00696134627379** |
| Probability mass with individual budget errors above `1e-5` | **0.0468806276181 (4.6881%)** |
| Probability mass with nonpositive consumption | 0 |
| Probability mass with original floor consumption `1e-5` | 0 |
| Probability mass violating borrowing minimum in saved fine savings policy | 0 |
| Maximum evaluated first-period Euler/KKT discrepancy | **0.000839922211416** |
| Euler/KKT-undefined positive probability mass | 0 |
| Positive probability mass with Euler/KKT discrepancy exceeding `1e-4` | 0.00831556130075 |
| Maximum absolute realized home-goods residual, periods 1–30 | **0.0043352276484** |
| Maximum source asset residual, periods 1–30 | 1.3536933352e-8 |
| Maximum source labor residual, periods 1–30 | 5.39561950674e-10 |

The table records unscaled model units unless the metric explicitly says otherwise. Results were produced by independent source-path evaluation in CI; this is **not** a claim of independent MATLAB recalculation or any nonlinear transition solve.

### Household feasibility: original un-interpolated states

In `backsolve_egm.m`, source `G0_shock` records the inherited 769×7 distribution, and `In_dis_exp` records the redistributed first-period fine-grid distribution. For an original wealth-grid node `a` with state `s`, the source portfolio maps

```text
k(a) = max(ell * a + k_const_hat*K0, 0)
b(a) = a - k(a)
ell = (1-k_const_hat)*q0*K0 / sum_{a>0,s} G0_shock(a,s)*a
a_post(a) = k(a) + b(a)*(x_pre/x1)
```

The independently checked *necessary* strict affordability inequality for positive consumption and next-period borrowing `a' >= bmin` is

```text
(1+r1)*a_post(a) + wage_income(s) + Profits1 - bmin > 0.
```

At the original MATLAB depreciation path there are **zero violating original states with positive mass**. The minimum equals `0.243344083414`. It is not appropriate to transplant the earlier Julia high-branch finding of 6.57% infeasibility to the original saved depreciation path. The latter features different macro income and price paths.

### Economic equations and where they fail

The original `solve_transition_bis.m` computes the scalar aggregate consumption from its asset-budget expression and determines the terms of trade by `fsolve(aggr,0.8)`. That scalar equation is **not** an external component of the 399-dimensional outer residual. With original saved `K,L,r,w,Profits,q,B`, at its returned `x_1`, the scalar difference is `+0.00211041927838`, not zero.

Furthermore, source `backsolve_egm.m` produces a revalued fine-grid distribution and fine consumption policy; the actual aggregate household home-goods demand exceeds the source aggregate-budget formula by `0.00222480837002`. Consequently, the *realized* goods residual is `−0.0043352276484`. The two distinct residual objects are not interchangeable. In levels:

```text
C_budget - C_goods_supply = +0.00211041927838
C_households - C_budget = +0.00222480837002
C_goods_supply - C_households = -0.00433522764840.
```

### Analytical exclusion of any scalar goods root at fixed authors' aggregate K/L

For the saved source calibration `theta=1`, the home-consumption expenditure factor `m=1+omega/(1-omega)` does not vary with the terms of trade. At fixed authors' aggregate capital, labor, wage, interest and profits, the original scalar residual has the exact functional form

```text
f(x) = C + E/x + D/x^theta_star,
E = (1+r_1)*BG_shock*x_pre/m < 0, D = d_s > 0.
```

With `theta_star > 1`, the *only* stationary point on positive prices is

```text
x_min = (theta_star*D/(-E))^(1/(theta_star-1)).
```

The function decreases before `x_min` and increases after it, so this is its global minimum, independently of any starting guess or numerical root scanner. For the actual saved MATLAB leveraged flexible-FX macro path, `x_min = 0.909875606122` and `f(x_min)=+0.00211041785147`. **There is no positive root of the source scalar goods equation at this fixed aggregate allocation.** MATLAB's returned `x_1=0.909858511411` lies near this residual minimum. Calling `fsolve` with a different starting guess cannot make the current fixed-K/L allocation satisfy this scalar equation. This does not rule out a different complete GE with different aggregates.

The original gross-portfolio shortfall (`−0.03188894`) is independently measured on the actual inherited MATLAB distribution; it differs slightly from the previously measured Julia-distribution shortfall (`−0.03189366`). The two cannot be conflated. Saving / aggregate total-wealth clearing to `1e-8` does not by itself prove that *gross* capital and foreign liabilities have their assumed aggregate levels.

At positive-mass saved fine-grid states, first-period budget residuals reach `0.00696135`, even though the stored consumption policies are positive and next-period savings do not breach the borrowing floor. Therefore neither initial household insolvency nor an active literal `1e-5` consumption floor explains the entire first-period discrepancy. A 4.6881% probability mass has first-period individual budget discrepancies exceeding `1e-5`. The independent audit exports every such positive-mass fine-grid state in `first_period_occupied_budget_errors.csv`, including its asset-grid index, productivity index, probability, saved consumption and savings, available resources, and signed budget residual. Pinpointing interpolation-versus-budget source mechanics requires comparing those exact states with the source constrained-savings and PCHIP expressions.

### Equilibrium status and what follows logically

1. **Source saved depreciation path: not a verified full equilibrium.** The actual first-period goods residual is nonzero and positive-mass budget residuals are non-negligible. The saved 399-equation convergence is insufficient.
2. **Original household affordability: passes the necessary pre-interpolation test.** At the authors' own low-price macro allocation, there is no positive-mass household with nonpositive available resources.
3. **Original first-period consumption floor: not active on saved positive-mass fine-grid nodes.** This does not rule out other constrained-policy interpolation problems.
4. **One major balance-sheet accounting inconsistency persists in the literal original portfolio mapping.** Correcting it would change the economic model relative to the literal authors' code and requires a separate GE recalculation. The already completed strict Julia GE still appreciates; it does not settle whether a different depreciation branch exists.
5. **No unrestricted conclusion about existence or uniqueness of a depreciation-side full GE.** The next computation, if any, must be justified after source-policy and market-closure reconciliation. A numerical solver failure cannot establish nonexistence.

### Limitations and next admissible step

The Stage 6 audit independently reconstructs the original saved shock-date affordability and first-period original policies, asset, labor, source scalar and realized home-goods residuals. It also reconstructs first-30-period realized home-goods residuals from the saved authors' aggregated consumption series. It does **not** yet certify all later-period household budget / Euler inequalities, fixed-FX price-setting conditions, the full terminal continuation, the asymptotic infinite-horizon equilibrium, or the economic cause of the discrepancy.

Before any 599-dimensional continuation, identify the positive-mass states with the largest household budget errors; compare the source fine-grid interpolation and constrained savings rule against the literal budget, and separately verify the analytic stationary minimum of the source scalar goods equation (rather than treating an `fsolve` exit status as a root). The new script reports the stationary candidate and its residual when the calibrated equations admit that characterization.

### Reproduction

From `FMR-2020/`:

```powershell
python -m pip install numpy scipy
python -m unittest discover -s test -p 'test_stage6_author_depreciation.py' -v
python scripts/stage6_author_depreciation_audit.py
```

The saved `.mat` in `results/author_reference/` is used unchanged. For the original MATLAB repository on Windows, a separate comparison can read the original file `C:\github\HANKSOME\exp_flex_lev16.mat` with `--source`, but is not required for CI. No MATLAB or Julia solver is invoked.
