# Stage 5 — Independent equilibrium-accounting and FX-root diagnosis

**Status (3 November 2026): definitive inherited gross-portfolio accounting error; conditional FX-price branching confirmed; full depreciation-side general equilibrium not yet found or excluded.**

Paper: de Ferra, Mitman, and Romei (2020), *Journal of International Economics* 124, 103303, [DOI](https://doi.org/10.1016/j.jinteco.2020.103303). Source: [HANKSOME MATLAB](https://github.com/kurtmitman/HANKSOME).

This document uses the existing **committed** Euler-polished 200-period Julia market roots. It does NOT infer paper equivalence from small market residuals, and it does not modify those archived roots.

## 1. Equation-level comparison (article vs authors' MATLAB vs Julia)

| Economic condition | Paper | Original authors' MATLAB | Current Julia | Assessment |
|---|---|---|---|---|
| Home/foreign relative price | (19)–(20), \`x=P_H/e\` with foreign foreign-price normalized to 1 | \`solve_transition_bis.m\`, \`t_guess\`; \`Figures.m\` nominal rate \`1/tot\` | \`w.terms\`, \`nominal_exchange_rate=terms_pre/terms\` (relative normalization) | Conventions consistent; **higher x means nominal appreciation** in floating-price regime |
| Foreign demand | (21), \`C_H^*=D*x^(-theta_star)\` | \`solve_transition_bis.m\`, lines ~139–144 | \`flexible_full_goods_equilibrium.jl\` goods residual, lines ~60–64 | Demand sign consistent |
| Home-goods resource | (22), \`Y=C_H+C_H^*+I+Psi\` | \`solve_transition_bis.m\`, \`net_cons\`; source's terminal is approximated | \`stage3_flexible_full_market_residual!\` (actual weighted household C) | Interior full home-goods clearing correct for the *implemented* prices |
| Net foreign credit-supply reversal | §5–6, unanticipated contraction of access to credit (B becomes less negative) | \`Main.m\` generates \`shock_b_agg\`; \`MIT_transition.m\` lines ~127–148 use \`new_shock_b_agg = Bcal + (Bpre-Bcal)*linspace(1,0,T)^4\` | \`credit_reversal_path\`, \`stage3_reversal_context\` | **Correct sign and t=41 timing**; unexpected B at t1 equals last prior B, but differs from previously expected B_t1 and all future paths |
| Initial portfolios | Household asset decomposition before unexpected currency movement | \`Transition/backsolve_egm.m\` lines ~40–55: \`k_i=max(lev*a_i+k_const,0)\` | \`revalue_foreign_debt!\` lines ~225–266 | **Source and legacy Julia both fail aggregate gross-capital clearing when k_const_hat>0** |
| Initial aggregate household budget | Household budget (6), with current returns and inherited wealth in home-good units | \`solve_transition_bis.m\`, lines ~129–135, uses \`q0*K+Bpre*xpre/x\` | \`stage3_fixed_prices!\`, lines ~77–87 | Budget *assumes* gross capital \`q0*K\`, contradicting the legacy actual individual portfolio allocation |
| Exchange-rate revaluation | §3 partial-equilibrium, household wealth effect of currency-price change | \`backsolve_egm.m\`, line ~127, \`Bgrid_fine/(t_guess(1)/tot_initial)\` | \`bonds*xpre/x\` | Correct **conditional** sign; not the source of appreciation |
| Union wage optimality | Paper (5), g'(L) factor | Source \`lab_prime\` omits \`1/varphi\` power though Main.m calibrates \`varphi=0.5\` | Julia includes nonlinear labor exponent | Source-level inconsistency already corrected in Julia, but this changes the source computational experiment |
| Marginal-cost/price setting | (14)–(15) and Rotemberg NKPC | \`solve_transition_fixed.m\` backward nonlinear marginal cost | \`stage3_fixed_prices!\` | No independent error established in existing pricing recursion |
| Capital/bond no arbitrage | (8)–(10), shared asset return under perfect foresight | Source returns from rental capital and q evolution | Julia computes same gross returns | No independent exogenous foreign interest-rate path to test in these data; do not invent one |
| Finite-horizon closure | Infinite-horizon model | Source imposes terminal policy and copies last consumption | Julia pins last-period L and x, not final K | **Approximation**, not a theorem; pre-shock 200/300 early paths stable, but surprise-path terminal tests remain open |

## 2. Exact aggregate portfolio inconsistency

The source/legacy Julia maps net household wealth \`a_i\` into gross capital and bonds:

\[
k_i^{raw}=\max\{\ell a_i+\widehat k K_0,0\},\qquad
b_i^{raw}=a_i-k_i^{raw},
\]
\[
\ell=(1-\widehat k)\frac{q_{-1}K_0}
{\sum_{i:a_i>0}\mu_i a_i}.
\]

But the GE budget block assumes the pre-shock aggregate gross holdings satisfy

\[
\sum_i\mu_i k_i=q_{-1}K_0,\qquad
\sum_i\mu_i b_i=B_{-1},\qquad
\sum_i\mu_i a_i=q_{-1}K_0+B_{-1}.
\]

For the actual 200-quarter saved Julia shock-date distribution, independently recovered by one re-evaluation of the archived corrected boom root, the measured gaps are:

| Experiment | Sum gross capital minus q-pre K | Sum bonds minus B-pre | Actual minus assumed home consumption at date 1 |
|---|---:|---:|---:|
| Flexible, capital intercept 0 | \(-4.62\times10^{-14}\) | \(-1.96\times10^{-7}\) | \(-1.46\times10^{-7}\) |
| Flexible, capital intercept 1/16 | **\(-0.0318936557\)** | **\(+0.0318934600\)** | **\(-0.00229820925\)** |
| Fixed, capital intercept 1/16 | **\(-0.0318936557\)** | **\(+0.0318934600\)** | **\(+0.00028322378\)** |

The accounting wedge can be derived without a model solve. Let \(\lambda_x=x_{-1}/x_1\), let \(\widetilde K=\sum\mu k^{raw}\), let \(\widetilde A=\sum\mu a\), and let \(R_1=1+r_1\). The true revalued total household wealth is

\[
A_1^{true}=\widetilde K+(\widetilde A-\widetilde K)\lambda_x,
\]

while the aggregate budget uses

\[
A_1^{agg}=q_{-1}K_0+B_{-1}\lambda_x.
\]

Therefore

\[
A_1^{true}-A_1^{agg}
=(\widetilde K-q_{-1}K_0)(1-\lambda_x)
+(\widetilde A-q_{-1}K_0-B_{-1})\lambda_x,
\]

and, with \(\chi=0.6\),

\[
C_{H,1}^{actual}-C_{H,1}^{budget}
=\chi R_1\big(A_1^{true}-A_1^{agg}\big).
\]

The numerical discrepancies **satisfy this independent wealth identity** within existing market-root error. This is an *actual model-accounting failure*, not the user-PC's ordinary numerical tolerance or a conjecture about code.

### Economic correction (explicitly opt-in, so past results remain comparable)

\`revalue_foreign_debt!(...;strict_portfolio=true)\` now multiplies the nonnegative gross capital portfolio by

\[
s_K=\frac{q_{-1}K_0}{\sum_i\mu_i k_i^{raw}},
\quad k_i^{strict}=s_Kk_i^{raw},\quad b_i^{strict}=a_i-k_i^{strict}.
\]

It then checks inherited aggregate wealth, gross capital, bonds, and mass. The common scaling preserves nonnegative capital and the relative pattern of gross exposure, but slightly changes the *exact* meaning of the intercept \(\widehat k\). This is a disclosed economic modeling choice, not an exact replication of the authors' original portfolio mapping. **The strict system must be solved again before claiming a new GE**, so existing four archived roots continue to use the legacy mapping by default.

No strict-portfolio 599-dimensional GE roots have yet been produced. The new strict-capital mapping is available as an opt-in and unit-tested component. Never overwrite the old outputs.

## 3. Two positive first-period *conditional* trade roots

With the existing aggregate K/L/q/w/r price paths held fixed, the t=1 aggregate trade equation is

\[
f(x)=C+\frac{E}{x}-\frac{D}{x^{\theta^*}}=0,
\quad
E=-\chi(1+r_1)B_{-1}x_{-1}>0,
\]

where \(\theta^*=3\), \(D\) is the foreign home-goods demand shifter, and \(C\) collects net home-good resources and capital, wages and profits excluding the inherited foreign bond's exchange-rate term.

The realized coefficient sets have **two strictly positive scalar roots**:

| Conditional-price experiment | Low root (depreciation) | High root (appreciation) | Low root preserves positive-mass household feasibility? |
|---|---:|---:|---|
| Flexible, intercept 0 | 0.68856213 | 1.17218198 | **Yes** |
| Flexible, intercept 1/16 | 0.69754172 | 1.15672207 | **No** |
| Fixed, intercept 1/16 | 0.79966750 | 1.01233838 | **No** |

These numbers were computed and tested in GitHub Actions #326 using the **saved actual Julia roots** and the actual household distribution. The conditional high root is close to the actual flexible-price solution; their difference under leverage is itself consistent with the aggregate portfolio discrepancy identified above.

At the saved K/L path, a source-style lower-price root *is not financially feasible* under the leveraged household model: some occupied low-wealth states fail to fund positive consumption and meet the borrowing limit. The source MATLAB uses a numerical consumption floor; treating that as an admissible household budget would be economically invalid.

**Important limitation:** Multiple scalar roots conditional on K/L do NOT prove multiple equilibria for the full 599-dimensional household/firm/asset/goods markets. The existence or nonexistence of a feasible depreciation-side full GE is not established.

The full 599-equation equations can only be re-solved on a depreciation-initialized branch after correcting the aggregate portfolio accounting and checking the shock-date feasible set. Do NOT simply flip the sign of exchange-rate conversion or assert multiple equilibria.

## 3.1. Decisive direct-state test: economic infeasibility versus grid interpolation

The initial Stage 5 low-root feasibility test evaluates **post-interpolation**
positive-mass fine-grid nodes, which is not by itself proof that the
*underlying* shock-revalued household states are insolvent.

To discriminate, commit `bc101fdac5f6d073218332784e4cda881d96413d`
adds `stage5_exact_original_state_support`, which calculates directly
for every positive-mass original state `(a_i,s_i)`:

```
k_i = capital_normalization * max(lev * a_i + kfixed, 0)
b_i = a_i - k_i
a_i_post = k_i + b_i * (x_pre / x_candidate)
available_i = (1+r_1)*a_i_post + w_1*L_1*s_i + profits_1 - borrowing_limit
```

**No interpolation of `a_i_post` precedes this check.**
GitHub Actions **#356 SUCCESS**, on saved 200-period roots, gives:

| Low-price conditional branch | Minimum disposable income after grid interpolation | Minimum income at original household state (BEFORE interpolation) | Mass of original states with negative resources |
|---|---:|---:|---:|
| Flexible, no capital intercept | +0.50202019 | **+0.50213965** | 0% |
| Flexible, capital intercept 1/16 | -0.46416575 | **-0.43030395** | **6.57173199%** |
| Flexible, capital intercept 1/16, capital-conserving correction | -0.46416575 | **-0.43128630** | **6.57173199%** |
| Fixed, capital intercept 1/16 (old fixed-price block) | -0.00092531 | **+0.02897273** | **0%** |

**Conclusion for the flexible-FX leveraged low scalar root conditional
on archived K/L/r/w/profits:** it violates the **model's own necessary
household budget/borrowing feasibility condition for 6.57% of original
positive-mass households**. This is not introduced by grid interpolation,
nor does it depend on EGM versus fsolve. The original source's
`max(c,1e-5)` does not repair a negative available budget.

**Qualification 1:** This result does **NOT** prove that every
depreciation-side full general-equilibrium path is impossible: another
K/L/price/return path changes household income and may be feasible.
A negative disposable resource for the same fixed aggregate path is
a conditional impossibility, not a global equilibrium nonexistence
theorem.

**Qualification 2:** The "fixed-FX lower scalar root is infeasible"
claim from the *interpolated* grid check must NOT be treated as an
economic theorem. Its original micro states actually have
**positive** disposable resources at the *archived* fixed-price block.
Furthermore, changing x in a **fixed** regime changes inflation,
Rotemberg marginal costs, factor prices and household returns. The
earlier scalar low root holding these prices fixed is not a valid
fixed-regime price-setting counterfactual.

## 4. Next numerical discrimination

1. Run the saved-root first-period audit with **strict_portfolio=true** and separately record gross capital/bond clearing, first-period household budget aggregation and lower-price household affordability. This is one household distribution/market residual evaluation, **not** another large Newton solve.
2. If the corrected portfolio removes the budget wedge but the low price remains infeasible, then the legacy portfolio-accounting discrepancy is real but **not the sole cause** of the depreciation failure.
3. Only then test controlled continuation in the full GE: preserve pre-shock wealth distributions and exogenous foreign credit path, warm-start the published-style low-price branch while enforcing feasibility (no consumption floor). Track whether it reaches a distinct valid root or fails.
4. Only if a price-response-consistent, budget-feasible equilibrium exists should the main paper conclusions (fixed/flexible consumption and welfare) be claimed as replicated.

Research commands (from \`FMR-2020/\`):

\`\`\`powershell
julia --project=. scripts/stage5_initial_fx_audit.jl
julia --project=. test/test_stage5_first_period_audit.jl
\`\`\`

No new 599/899-dimensional root solve is required by these commands. CI snapshots of the exact numeric diagnostics are attached as the Stage 5 first-period accounting artifact.

## 5. Strict-portfolio one-pass counterfactual: realized test results

GitHub Actions run #345 and #346 initially FAILED because a test
incorrectly treated a corrected initial portfolio as a new equilibrium
**without re-solving household consumption, prices, and the asset and
goods markets**. This is the wrong numerical claim: the two portfolio
accounting constraints can clear exactly while the old aggregate path
ceases to be a GE root. The test has been corrected accordingly.

At the archived (unchanged) Euler-polished GE aggregate price and
capital paths, enabling `strict_portfolio=true` gives:

| Diagnostic | Flexible intercept 1/16 | Fixed intercept 1/16 |
|---|---:|---:|
| Gross capital minus inherited qK | approximately 3.2e-14 | approximately 3.2e-14 |
| Gross foreign bonds minus inherited B | within historical 2e-7 asset-market error | within historical 2e-7 |
| Weighted household C_H minus old aggregate-budget C_H | **-0.00223858** | **+0.000273884** |
| Full 599-system scaled market residual at OLD price path | **1.17791e-4** | **1.54654e-5** |
| Low conditional price root feasible at OLD macro K/L | **NO** (minimum disposable resources **-0.46417**) | **NO** (minimum **-0.0009253**) |

Thus the corrected gross portfolio **does not automatically restore**
the paper's depreciation branch, nor does it give a new GE solution.
In particular, the lower scalar root is still not affordable at the
old aggregate K/L path. A different endogenously adjusted K/L path
might alter household budgets; the low-root full GE is therefore
**unproven, not impossible**.

The original positive-mass household feasibility checks and the
conditional scalar price-root values remain unchanged as benchmarks.
The old roots must NOT be silently renamed "strict-portfolio
equilibria". A strict-portfolio 599-equation solve would be an
entirely distinct research computation and should be undertaken
only if needed to distinguish the remaining equilibrium-selection
hypothesis.

## 6. Controlled experiment: one strict-portfolio 599-market solve

After establishing the exact inherited-capital accounting failure and
that simple reallocation at old macro prices **does not** clear the
599-equation GE, the next causally informative full simulation is
**only the leveraged flexible-FX sudden stop**, holding the corrected
Euler-polished credit-expansion path, foreign-credit shock, and
predetermined household wealth distribution fixed.

A dedicated opt-in mode is available in the existing runner:

```powershell
Set-Location C:\github\22p8z0omwz9xu
git pull --ff-only
Set-Location FMR-2020
julia --project=. scripts\run_economic_replication.jl --stage5-solve-strict-flexible
```


It loads/rechecks archived corrected boom, reconstructs the new shock
context, warm-starts from the archived leveraged-flexible
Euler-polished 599-dimensional path, and re-solves the SAME market and
Euler/KKT equations with `strict_portfolio=true` inside every
household distribution and finite-difference Jacobian evaluation.
The new result (only if converged) is saved as:

`economic_flexible_leveraged_euler_polished_STRICT_PORTFOLIO_MARKET_ROOT_UNCERTIFIED.csv`.

This separate CSV never replaces the legacy portfolio roots. The
runner explicitly reports impact x, implied nominal appreciation/
depreciation, market residual, Euler error, and terminal capital gap.

Interpretation of results:
- If the new strict equilibrium **still appreciates**, this specific
  aggregate portfolio accounting error does not explain the published
  exchange-rate sign. The low-branch GE existence/selection and model
  assumptions remain the central issue.
- If the strict equilibrium **depreciates** and satisfies household
  budgets and all markets, then the gross portfolio accounting
  inconsistency is sufficient to switch the GE response in this
  calibration. Only then evaluate the new fixed-FX equilibrium and
  distributional welfare counterfactual.
- If the solver does not converge, the experiment is **inconclusive**;
  do not interpret nonconvergence as proof of no equilibrium branch.
  Preserve the exception/trace for targeted diagnosis.

This justified one-economy GE experiment **has now been executed on
the user's Windows PC**. Its outcome is documented below.

## 7. Strict-portfolio flexible-FX root: actual user-PC numerical result

The user ran from `FMR-2020/`:

```powershell
julia --project=. scripts\run_economic_replication.jl --stage5-solve-strict-flexible
```

The archived 200-quarter Euler-polished credit expansion was independently
rechecked (scaled max residual **2.558332579e-8**). The single
leveraged flexible-FX sudden-stop root was recomputed with inherited
gross capital and foreign bonds clearing exactly, not just at an
unchanged old aggregate path. The 599-equation solver required 3
iterations and 605 residual evaluations.

| Diagnostic | Archived leveraged flexible root (legacy portfolio) | Recomputed strict portfolio |
|---|---:|---:|
| Terms before the shock | 1.0271423614 | 1.0271423614 |
| Impact terms x=P_H/P_F | **1.1584047924** | **1.1561011378** |
| Nominal depreciation, 100 log(x_pre/x_1) | **-12.0263 log points** | **-11.82727 log points** |
| Maximum scaled market residual | 2.650952e-7 | 1.230285e-7 |
| Maximum raw asset-market discrepancy | 7.845209e-6 | 3.640897e-6 |
| Maximum raw goods-market discrepancy | 1.520263e-7 | 6.674653e-8 |
| Maximum household budget error, transition | 4.44e-16 | 4.44e-16 |
| Maximum Euler/KKT violation | 9.997379e-6 | 9.997586e-6 |
| Distribution mass conservation error | 1.25e-14 | 1.22e-14 |
| Terminal capital level gap | 0.010511 | 0.010456 |
| Terminal household budget max | 0.011091 | 0.011342 |

The root remains a **substantial NOMINAL APPRECIATION** and thus
benefits short foreign-currency debt positions under the current
revaluation convention. Fixing the inherited household gross
capital/bond accounting by itself changes the price response by only
about **0.20 log points** and **DOES NOT** reproduce the
depreciation-led balance-sheet mechanism in the paper.

This result is **from the user's logged local execution**, not a
full-GE root independently rerun on GitHub Actions.
The separately saved root is at:

`results/economic_reproduction/economic_flexible_leveraged_euler_polished_STRICT_PORTFOLIO_MARKET_ROOT_UNCERTIFIED.csv`

The result has not passed all terminal, asset-grid or
200-vs-300 sudden-stop horizon checks, and the existence of a
**different full depreciation-side GE** remains open. The fact
that the previously computed *conditional* low trade-price root
makes 6.57% of inherited original-mass households insolvent
under HIGH-BRANCH K/L/r/w/profits establishes only
**conditional** impossibility, not a global nonexistence result.

A separate resumption bug was corrected in
`scripts/run_economic_replication.jl`: when rechecking this
strict-portfolio CSV, flexible-FX residuals must receive
`strict_portfolio=true`, just as fresh residual/Jacobian
calls do. It does not change the first execution's valid
strict-map solve, but prevents a repeat run from testing the
wrong household portfolio mapping.

**Decision:** Do not repeat the already converged
strict-portfolio root. Investigate the original MATLAB's
lower branch at its OWN prices, wages, and capital/labor paths,
using the continuous pre-interpolation household budget test.
Only if that is a feasible all-market root should a
depreciation-side full-GE continuation be attempted.
Any second-root search must vary the endogenous macro paths,
rather than substitute the low scalar price root while
holding the appreciation-branch prices/incomes fixed.

## 8. Unresolved limitations

- Which full GE branch, if any, corresponds to the original paper's depreciation scenario? **Not established.**
- Is the authors' archived low-price MATLAB approximation feasible after removing consumption floors and enforcing all goods and asset equations? **Not established in the new strict system.**
- Does exact initial gross-capital portfolio conservation by itself restore the depreciation-led FX mechanism? **No. The user-PC strict GE still appreciates by 11.83 log points.**
- Distributional welfare, the price-policy counterfactual, grid convergence and sudden-stop horizon convergence must be evaluated **after** the economic mechanism is resolved.
