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

## 5. Unresolved limitations

- Which full GE branch, if any, corresponds to the original paper's depreciation scenario? **Not established.**
- Is the authors' archived low-price MATLAB approximation feasible after removing consumption floors and enforcing all goods and asset equations? **Not established in the new strict system.**
- Could a correctly constructed pre-shock portfolio that strictly clears gross capital restore the original paper's FX mechanism? **Needs numerical test; the portfolio wedge is small relative to the ~12-log-point appreciation.**
- Distributional welfare, the price-policy counterfactual, grid convergence and sudden-stop horizon convergence must be evaluated **after** the economic mechanism is resolved.
