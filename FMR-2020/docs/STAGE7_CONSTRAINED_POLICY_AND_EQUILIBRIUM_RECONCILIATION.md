# Stage 7 — Constrained Policy Root Cause and General-Equilibrium Reconciliation

**Status: source household budget defect identified and its impact measured; no new depreciation-side complete GE certified.**

Research object: Sergio de Ferra, Kurt Mitman and Federica Romei (2020), *Journal of International Economics* 124, 103303, [published article](https://doi.org/10.1016/j.jinteco.2020.103303). The household budget is the paper's equation (1), expressed in the source code's home-good numeraire and two-good expenditure convention. The actual saved, leveraged floating-exchange-rate, shock-date experiment is **\`results/author_reference/exp_flex_lev16.mat\`**. All numbers below use that original saved MATLAB workspace; none are inferred from the Julia appreciation-branch root.

## 1. Source-level root cause (proven, not hypothetical)

Original MATLAB: \`HANKSOME/Transition/backsolve_egm.m\`.

- Lines **68–71** construct **coarse-grid** constrained consumption \`c_constrained = max(((1+r_t)*Agrid + labor_income + profits - bmin)/m, 1e-5)\`.
- Lines **100–101** generate fine-grid constrained consumption by \`interp1(Agrid,c_constrained(:,j),Agrid_fine,'pchip')\`.
- Lines **103–107** separately assign next-period savings \`dec_temp_fine=bmin\` whenever the constraint binds.
- Lines **126–152** revalue inherited gross foreign positions and interpolate the original distribution onto the fine grid.
- Lines **160–166** aggregate **actual policy consumption**, which need not equal the consumption used by the nested scalar trade-price solver.
- \`Transition/solve_transition_bis.m\` flexible-price lines **128–158** determine trade-price \`x_t\` with \`fsolve(aggr,0.8)\` and compute separate aggregate-budget and goods-resource consumption.

The constrained-household budget, using exact source units/timing, is

\`\`\`text
m(x_t) c_t(a,s) + a'_t(a,s)
    = (1+r_t) a + markup*(1-alpha)*Y_t*ex(s)/ex_mean + Profits_t,
m(x_t)=1+x_t^(theta-1)*omega/(1-omega),
a'_t(a,s)=bmin at the binding constraint.
\`\`\`

Hence **correct fine-node constrained consumption** is affine in current wealth:

\`\`\`text
c_exact(a,s) = ((1+r_t)*a+income_t(s)-bmin)/m(x_t).
\`\`\`

A PCHIP of the *unfloored affine coarse-grid budget* reproduces this to numerical precision, but a PCHIP of \`max(c_exact,1e-5)\` generally does not. The source's cubic coarse/fine asset grids include a negative below-borrowing-limit tail, so the floor is active at **178 original coarse asset–productivity nodes** across seven states. It changes neighboring PCHIP derivatives even where the floor is not active at the fine-grid consumption node being evaluated.

**Exact source-data discrimination:**

| Diagnostic on the original shock-date positive-mass binding states | Verified maximum / value |
|---|---:|
| Original saved fine consumption minus reconstructed source PCHIP of floored coarse values | **2.22044604925e-16** |
| PCHIP of *unfloored* affine coarse budget minus exact fine-node budget | **1.110223024625e-16** |
| PCHIP with coarse floor minus PCHIP without floor (absolute maximum) | **0.004176807764276** |
| Weighted PCHIP-floor consumption difference on binding states | **−0.00009567180543369** |
| Original positive-mass states with budget discrepancy above \`1e-5\` | **4.68806276181%** |
| Positive-mass nonbinding states with such discrepancy | **0** |

Thus **the coarse-grid floor followed by PCHIP is the source of the entire substantive first-date individual budget violation**; the standalone PCHIP operator, source income timing, and the unfloored affine formula pass the targeted algebraic tests. This conclusion is grounded in the saved MATLAB policies, not a guess about EGM implementation.

The worst source fine-grid state is **(asset index 116, productivity index 1)**, with \`a=-0.2487936594115\` and positive mass \`0.00762231849589\`. Source consumption is \`0.1351266584392\`, while the exact binding-budget value is \`0.1393034662034\`; its signed budget residual is \`−0.006961346273794\`. In that productivity state, 28 coarse asset nodes have floor-active constrained consumption. The worst fine node lies just above the cutoff, so adjacent floored coarse knots distort PCHIP slopes even though the two enclosing coarse consumptions are themselves unfloored.

This is **not** the same as financial infeasibility of the original pre-interpolation households. Stage 6 found no positive-mass state with nonpositive available resources at the original authors' own prices and wages, and this finding remains intact. A shock-revalued current wealth \`a\` can lie below the *next-period borrowing limit* \`bmin\` even when consumption is strictly feasible. The coarse negative grid tail was retained; we did not trim, coarsen, or change the asset grid.

## 2. Counterfactual correct binding budget, original macro path fixed

The new independent \`scripts/stage7_constrained_policy_audit.py\` keeps the original author-saved \`K,L,x,r,w,profits\`, distribution \`mu\`, savings \`a'\`, foreign-credit path, and portfolio mapping **exactly as stored**. ONLY positive-mass fine-grid cells with \`a'=bmin\` receive

\`\`\`text
c_corrected = ((1+r)*a + income(s) - bmin)/m.
\`\`\`

It explicitly rejects nonpositive exact consumption and does not introduce any numerical floor. Nonbinding policies are unchanged. It reconstructs source PCHIP with SciPy and requires machine-level agreement. **This is a partial-equilibrium policy-identity counterfactual; not a new policy optimum or general equilibrium.**

| Diagnostic, date 1 | Original saved MATLAB | Exact binding-budget counterfactual |
|---|---:|---:|
| Maximum occupied individual budget error | **0.006961346274** | **5.830447236e-12** |
| Aggregate household home-good consumption | 1.336082854225 | **1.336178526031** |
| Goods supply available for home consumption | 1.331747626577 | unchanged |
| Realized goods-market residual, supply less consumption | **−0.004335227648398** | **−0.004430899453832** |
| Maximum measured first-date Euler/KKT discrepancy | 0.000839922211416 | 0.000839922211416 |
| Minimum positive-mass resources before source interpolation | +0.243344083414 | unchanged |
| Source distribution mass | 1 | unchanged |

\`\`\`text
C_corrected - C_original = +0.00009567180543346.
Goods residual change   = −0.00009567180543346.
\`\`\`

**The budget correction fully repairs the detected micro-budget error but moves the original fixed-aggregate goods-market residual further away from zero.** Euler/KKT has NOT been reoptimized: the unchanged maximum is a diagnostic, not a certification of the corrected consumption/savings pair as a jointly optimal policy. No claim is made that the original solver's price \`x_1=0.909858511411\` clears the goods market after this change.

### All 199 interior-period original MATLAB policies

A second read-only script, \`scripts/stage7_full_transition_budget_audit.py\`, reconstructs the **entire original fine-grid distribution propagation for dates 1–199** from the saved household savings and the source Markov transition, with linear lottery allocation on neighboring grid nodes. It independently checks the reconstructed original household aggregate consumption against the original saved \`c_t_exp\` at every date. No macroeconomic root finding or fresh household Bellman solve occurs.

Observed CI results:

| Diagnostic | Value |
|---|---:|
| Interior dates audited | **199** |
| Dates with positive-mass budget residual above \`1e-5\` | **1 (the surprise-shock date)** |
| Maximum original individual budget residual across all dates | **0.006961346274** |
| Maximum after exact binding-budget correction across all dates | **6.960654275e-12** |
| Maximum affected mass in any date | **0.04688062761805** |
| Maximum forward-distribution mass error | **1.55431223447522e-15** |
| Largest single-date aggregate consumption change | **+0.00009567180543368** |

The original MATLAB next-period savings and lottery distribution are left unchanged. The correction is a budget-identity calculation, not a replacement equilibrium policy. No terminal-date \`t=200\` savings policy is stored in the original source, so the audit correctly stops at \`t=199\`.

## 3. Exact economic decomposition of the aggregate household consumption discrepancy

The original source has **two distinct goods-market failures**. With \`C_budget\` from \`solve_transition_bis.m\`, \`C_goods\` from the home-goods resource constraint, and \`C_households\` from the saved fine-grid policies:

\`\`\`text
C_budget - C_goods       = +0.00211041927838
C_households - C_budget = +0.00222480837002
C_goods - C_households  = −0.00433522764840.
\`\`\`

Let \`lambda=x_pre/x_1\`, \`R=1+r_1\`, \`m=m(x_1)\`, and define

\`\`\`text
A_post_model  = q0*K0 + B0*lambda,
A_post_micro  = sum(mu_post * a_grid),
A'_supply     = B_1 + q_1*K_2,
A'_micro      = sum(mu_post * savings),
E_B           = sum(mu_post*(m*c + savings - R*a - income)),
income_macro  = w_1*L_1 + Profits_1,
income_micro  = sum(mu_post*income(s)).
\`\`\`

Summing micro budgets **exactly** gives

\`\`\`text
C_households - C_budget
 = [R*(A_post_micro-A_post_model)
    +(income_micro-income_macro)
    -(A'_micro-A'_supply)+E_B] / m.
\`\`\`

Insert inherited gross-portfolio revaluation \`A_revalued=sum(mu_pre*(k+lambda*b))\`, and split

\`\`\`text
A_post_micro-A_post_model
 = (A_post_micro-A_revalued)+(A_revalued-A_post_model).
\`\`\`

Independent numerical contributions on the *actual original saved allocation* are:

| Mechanism | Contribution to \`C_households-C_budget\` |
|---|---:|
| Original gross capital/bond portfolio mismatch after foreign-liability revaluation | **+0.002320487927977** |
| Difference between post-shock wealth on grid and revalued original household wealth (distribution transport) | 0 to numerical precision |
| Labor income aggregation | +2.823075107017e-12 |
| First-period asset-market residual | −7.755345166061e-09 |
| Original constrained-household signed budget residual | **−0.00009567180543341** |
| **Sum** | **+0.002224808370021** |

Numerical decomposition error: **1.914701036609e-15**.

The source's literal portfolio formula allocates gross capital **0.0318889358579 below** required \`q_0 K_0\` and foreign bonds **0.03188894913191 above** required inherited \`B_0\`, while aggregate pre-shock *net* wealth still clears to \`1e-8\`. Depreciation makes the error economically relevant via \`A_revalued\`. This is a source micro–macro accounting discrepancy, not evidence that foreign-debt revaluation itself should change sign.

The exact constrained-consumption repair removes the **negative** individual-budget contribution, leaving \`C_corrected-C_budget=+0.002320480175453\`. Thus household budget repair **increases**, rather than eliminates, the aggregate consumption discrepancy. The gross-portfolio revaluation wedge is an order of magnitude larger than the constrained consumption correction.

These are additive accounting contributions **at the original saved aggregate allocation**. They must not be called additive effects of separate *general-equilibrium re-solves*; changing the original gross portfolio would change the shock-date wealth distribution, optimal household policies, and GE prices.

## 4. Independent scalar goods-market closure remains invalid

Stage 6 independently established, at the saved authors' aggregate \`K,L,r,w,profits\`:

\`\`\`text
f(x)=C+E/x+D/x^theta_star,
E<0, D>0, theta_star>1.
x_min=(theta_star*D/(-E))^(1/(theta_star-1))
     = 0.909875606122
f(x_min)=+0.00211041785147 > 0.
\`\`\`

Thus the original nested MATLAB \`fsolve(aggr,0.8)\` does **not** have a positive scalar zero at that fixed aggregate allocation. It returns a stationary point of the squared residual, not a goods-market root. No change in the MATLAB \`fsolve\` initial guess or tighter solver tolerance can repair this fixed-K/L failure. By itself this does **not** exclude a different *full endogenous K/L/x* depreciation equilibrium.

Moreover, source-level constrained consumption accounts for only \`0.000095672\` of the realized goods gap \`0.004335228\` in magnitude, and in the opposite direction from the repair needed to reduce the goods-market discrepancy. **It cannot be the source of the published exchange-rate mechanism's sign reversal in the converged Julia appreciation equilibrium.**

## 5. Corrected household solver: existing Julia implementation already enforces the budget

No replacement Julia EGM solver has been created because the existing relevant code is economically correct on this margin:

- \`src/sudden_stop_equilibrium.jl\`, function \`stage3_household_backward!\`, computes \`cons=(gross_old*a+income-borrowing_limit)/m\` **directly at the current fine asset** when binding (approximately lines 318–325).
- The nonbinding policy is interpolated on endogenous EGM asset knots and savings is computed from the *same budget identity*.
- The code rejects infeasible or unsupported positive-mass states, rather than silently replacing negative consumption by \`1e-5\`.
- \`src/transition_household.jl\` applies the same direct fine-budget principle.
- \`src/transition_distribution.jl\` and the Stage 3 forward-distribution block propagate positive-mass savings using the full fine grid; previous tests of these operators and saved Julia equilibrium roots are archived.

The Stage 7 independent Python counterfactual now verifies precisely the economic reason for the Julia design. Changing it again would be unnecessary and risk breaking previously passing Euler/market tests. A full MATLAB clone with the floor removed would still require separate proof that off-grid **zero-mass** nodes do not contaminate backward-EGM continuation; merely removing the numerical floor on the original full negative tail could create nonpositive consumption in unsupported states. We therefore **do not modify** the original MATLAB source or falsely describe the fixed-savings partial-equilibrium counterfactual as a new fully optimal household recursion.

## 6. Equilibrium decision: outcome B, with a limited existing appreciation comparison

**Outcome B: the original saved depreciation path fails full-equilibrium certification, but the existence of another depreciation-side full GE is unresolved.**

Evidence: original inner scalar equation has no fixed-allocation positive root, realized household goods do not clear, and gross portfolio accounting is inconsistent with aggregate inherited capital and bonds. The first-date household-budget bug is now fully explained and locally repaired; it is *not* sufficient to cure the market residual.

Existing Stage 5 Julia strict-portfolio 599-equation solution has \`x_1=1.1561011378\`, a nominal appreciation of approximately **11.83 log points**, with small finite-horizon market and Euler errors. This result is **already computed on the owner PC**, not recomputed by Stage 7. It remains subject to terminal and grid certification. Do not confuse an appreciation root under corrected Julia policy/accounting with proof that no depreciation-side equilibrium exists.

**Decision on another full nonlinear solve:** No indiscriminate 599/899-variable search is warranted. To establish a *distinct* depreciation equilibrium one must jointly vary original \`K,L,x\`; the fixed-K/L source scalar equation has no root, so neither scalar multistart nor an ex-post consumption patch can succeed. A controlled continuation may be informative only after confirming identical shock-date initial distributions, asset timing, portfolio convention, full household budget/Euler/KKT, all market equations, and terminal continuation. The existing \`scripts/evaluate_author_explicit_goods.jl\` can evaluate the original candidate with the current full Julia operator without solving; prior Stage 3 records report a nonzero full-599 residual at the literal original source candidate. Repeating the same audit alone is not progress.

Given that this Stage 7 run was executed in remote lightweight CI with **no local Julia or MATLAB runtime and no expensive full-GE authorization**, it does not claim new GE convergence, exhaustive multiple-root detection, or an infinite-horizon nonexistence theorem. The next *economically justified computational experiment*, if carried out on the owner PC, is a **single controlled low-price 599-dimensional continuation using existing solvers, with feasibility and all-market residuals recorded at every stage**. Separate strict-portfolio and source-portfolio systems must not be mistaken for multiple equilibria of one model.

## 7. Reproduction, source integrity, and tests

All runs read the tracked original author's \`.mat\` only and produce separately named CSV / Markdown artifacts:

\`\`\`powershell
Set-Location C:\github\22p8z0omwz9xu
git switch research/fmr-2020-stage3-sudden-stop-exchange-regimes
git pull --ff-only
Set-Location FMR-2020
python -m pip install numpy scipy
python -m unittest discover -s test -p 'test_stage7_constrained_policy.py' -v
python -m unittest discover -s test -p 'test_stage7_distribution_transport.py' -v
python scripts/stage7_constrained_policy_audit.py
python scripts/stage7_full_transition_budget_audit.py
\`\`\`

CI: \`.github/workflows/fmr-stage6-author-audit.yml\` also runs Stage 6 regression tests, \`py_compile\`, all new unit tests and both real original-source audits, archiving:

- \`results/stage7_constrained_policy/stage7_metrics.csv\`
- \`results/stage7_constrained_policy/stage7_binding_household_comparison.csv\`
- \`results/stage7_constrained_policy/stage7_full_transition_periods.csv\`
- \`results/stage7_constrained_policy/stage7_full_transition_summary.csv\`
- \`results/stage7_constrained_policy/stage7_audit.md\`

No original MATLAB workspaces, old Julia checkpoints or previously converged equilibria are overwritten.

**Central conclusion:** The first-period binding-household policy error is exactly traced to the *coarse consumption floor before PCHIP*; fixing it restores household budgets to roughly \`1e-11\` for all 199 interior dates, **but does not close the general-equilibrium goods market**. The remaining major components are separately and quantitatively established as a nonzero nested scalar goods equation and a gross-portfolio revaluation accounting wedge. Whether a distinct fully feasible depreciation branch exists is a controlled full-system equilibrium question, not an EGM interpolation question.
