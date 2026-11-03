# Stage 4 — Economic-direction audit: current GE roots do NOT reproduce the paper's main FX mechanism

Status: **Substantive replication discrepancy, unresolved.** This note is not a claim that the paper is wrong or that the Julia system has no equilibrium. It documents a mismatch between the currently converged market roots and the mechanism described in de Ferra, Mitman and Romei (2020), *Journal of International Economics* 124, 103303, Section 6, Figures 4–5.

## Evidence from actual committed 200-date Euler-polished Julia root files

Source files under `results/economic_reproduction/`:

- `economic_credit_expansion_euler_polished_MARKET_ROOT_UNCERTIFIED.csv`
- `economic_flexible_unleveraged_euler_polished_MARKET_ROOT_UNCERTIFIED.csv`
- `economic_flexible_leveraged_euler_polished_MARKET_ROOT_UNCERTIFIED.csv`
- `economic_fixed_leveraged_euler_polished_MARKET_ROOT_UNCERTIFIED.csv`
- `PROVISIONAL_economic_comparison.md`, `PROVISIONAL_figure4_5_economic_comparison.csv`

The current calibration is `chi=0.6`, `gamma=theta=1`. Under the implemented CES demand and price convention, `x=P_H/P_F`; under the flexible nominal exchange rate `e∝1/x`. The index of physical composite consumption is `C ∝ C_H^chi C_F^(1-chi)`. All quoted consumption comparisons are **100 times differences of logarithms, NOT CEV welfare**.

| Object | Continuing boom at calendar t=41 | Leveraged flexible sudden stop at t=1 | Leveraged fixed sudden stop at t=1 |
|---|---:|---:|---:|
| Net foreign asset supply B | -6.38201418 | -6.28100574 | -6.28100574 |
| Capital K | 29.59393015 | 29.59393015 | 29.59393015 |
| Labor L | 0.98600603 | 0.92987471 | 0.97389384 |
| Output Y | 3.02967599 | 2.91300509 | 3.00468990 |
| Terms x=P_H/P_F | 1.02614707 | 1.15840479 | 1.01200062 |
| Home consumption C_H | 1.46841413 | 1.54161636 | 1.42651982 |
| Foreign consumption C_F | 1.00453924 | 1.19054385 | 0.96242596 |
| Investment I | 0.62168518 | 0.71417011 | 0.58803384 |

The **last pre-shock realized** terms at boom date t=40 are `x_pre=1.0271423614420732`. The real-foreign-debt revaluation multiplier is `x_pre/x_1`; this is explicitly implemented as `bonds/(x_1/x_pre)` in `revalue_foreign_debt!`.

- Floating: `x_pre/x_1=0.8866869061`. Equivalently, the floating nominal-exchange-rate log depreciation is **-12.0263 log points**, i.e. a **12.03 log-point appreciation**. A negative foreign liability becomes **less negative**, not more onerous.
- Fixed: `x_pre/x_1=1.0149621853`. The real home-good value of foreign liabilities rises by approximately **1.50%** due to domestic deflation, despite the constant nominal exchange rate.
- Floating output relative to the continuing boom: **-3.9270 log points**. Fixed output: **-0.8281 log points**. Fixed output is **3.0989 log points above floating**, opposite the paper's description of a flexible-output boom and a fixed-output recession.
- On-impact composite consumption relative to continuing boom: floating **+9.7142 log points** versus fixed **-3.4498 log points**.
- On-impact fixed-minus-floating composite consumption: **-13.16395 log points**.
- On-impact leveraged-minus-unleveraged floating composite consumption: **-0.61136 log points**. This last comparison agrees in direction with the paper's idea that household foreign-currency leverage worsens the shock, but the first-order exchange-rate response reverses the model's intended balance-sheet channel.
- Over 30 quarters, the mean fixed-minus-floating composite consumption is **-1.25190 log points** and the ratio of 30-quarter sums of composite-consumption indexes fixed versus floating is **-1.27554%**. Neither is household welfare.

## Important counterfactual classification: "no leverage" is not no foreign debt

The existing `flexible_unleveraged_euler_polished` name is inherited from the numerical pipeline, **but it is an imprecise economic label**. The source call `stage3_reversal_context(...; leverage_hat=0.0)` leaves `foreign=true`. The underlying gross-asset allocation is

    lev = (1 - leverage_hat)*K_initial*q_initial / rich_wealth
    capital_i = max(lev * a_i + leverage_hat * K_initial, 0)
    foreign_bonds_i = a_i - capital_i.

Therefore `leverage_hat=0` only removes the FIXED CAPITAL INTERCEPT given to very low-wealth households. It does **not** eliminate foreign-currency liabilities or set foreign-currency borrowing to zero. When `lev>1`, even positive-wealth households hold net foreign-currency debt. The observed -0.61136 log-point contrast is a change in the **distribution of capital exposure/leverage**, not the effect of introducing foreign debt versus having no debt at all. The report preserves legacy CSV field names for compatibility, but explains the correct structural interpretation.

The current roots **do** satisfy the user-PC's full finite-horizon 599 market equations with scaled max residual <= 3.2410e-7 and Euler/KKT max < 1e-5. The 200-vs-300 boom's initial-40-period difference is <= 6.1931e-7. **These successes do not resolve this direction discrepancy**, and extending ONLY the initial boom horizon cannot certify the subsequent sudden-stop response.

## Published comparison and research implication

Section 6 of the 2020 article says the flexible FX response entails a deterioration in the terms of trade, foreign-currency debt revaluation, consumption contraction and an output boom. Fixed FX instead restrains terms-of-trade depreciation and debt revaluation, mutes the consumption loss and contracts output. It reports utilitarian welfare changes of -0.81% CEV under floating and -0.34% under fixing. Those published **welfare** statistics cannot be compared with Julia **aggregate consumption** changes directly.

The **currently committed Julia equilibria exhibit the opposite FX price movement, output ranking and aggregate consumption ranking**. At this stage do not claim reproduction of Figures 4/5, the paper's mechanism, or its welfare conclusions.

## Next numerical research decision

**Stop refining already-small market residuals and stop prolonging the boom.** Prioritize distinguishing **(A)** a mis-specified economic equation or timing and **(B)** selection of a different nonlinear equilibrium branch. Do not modify the code simply to force the published signs.

1. **Independent equation/timing audit**, anchored to paper Section 2 and original source code: verify that the net foreign asset path really constitutes a credit supply *contraction*; that the sign convention of x=P_H/P_F, the foreign export demand D*x^(-theta*), the foreign-currency liability rule `x_pre/x_1`, shock-date inherited assets, q_0, nominal exchange rate, and t=1 asset market agree with the paper. Revaluation's sign in existing Julia code is mechanically consistent given x; the cause of the appreciation must be diagnosed, not silently reversed.
2. **Check alternative roots only if equations/timing pass**. Use a controlled homotopy from the no-surprise credit path to the full announced reversal, and trace impact x, consumption and market residuals. Warm-start with previous converged roots, preserve all previous data, and test at least one separate depreciation-side initial guess. A genuine alternative requires independently converged market/budget/Euler/KKT conditions; a poor residual minimizer is NOT a root. If different valid roots exist, document an economically motivated equilibrium-selection criterion.
3. **Only after sign/equilibrium consistency**, compare fixed-versus-floating welfare and heterogeneous outcomes, target asset-grid refinement and sudden-stop 200-vs-300 sensitivity, and finish Figures 4/5.

The lightweight reporting script `scripts/report_polished_economic_results.jl` now generates `PROVISIONAL_mechanism_sign_audit.md` and explicitly reports any opposite signs; it deliberately does not equate convergence to successful replication.

## Reproducibility

From `FMR-2020/`:

```powershell
julia --project=. scripts/report_polished_economic_results.jl
```

No 599/899-dimensional nonlinear GE recomputation is needed to reproduce this audit; all required four-root CSVs are present in this branch.
