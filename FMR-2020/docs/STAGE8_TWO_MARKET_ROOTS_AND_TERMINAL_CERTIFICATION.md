# Stage 8 — Two strict-portfolio finite-horizon market roots: owner-PC evidence and terminal certification gap

**Status as of November 3, 2026: TWO DISTINCT 200-period strict-portfolio market-clearing candidates have been obtained from low- and high-price initialization. They are not yet certified infinite-horizon economic equilibria.**

Target: de Ferra, Mitman and Romei (2020), *Journal of International Economics* 124, 103303.
This report supplements and updates \`STAGE8_DEPRECIATION_BRANCH_EQUILIBRIUM.md\`. See the latter for the original-source initialization, detailed Stage 6–7 reconciliation, and actual lower-branch 599-Jacobian construction.

## 1. New evidence actually executed on the owner's Windows 11 / Julia 1.12.6 PC

The owner ran:

\`\`\`powershell
julia --project=. scripts\stage8_depreciation_continuation.jl --solve
julia --project=. scripts\stage8_depreciation_continuation.jl --resume
\`\`\`

Both commands used the same pre-shock Julia Euler-polished credit-expansion path (saved prior equilibrium), the same shock at period 41, strict inherited gross-capital/bond accounting, the same credit-reversal path, and the existing \`stage3_flexible_full_residual!\` with \`polish_households=true\`, \`strict_portfolio=true\`.

The first call accepted five full-system steps before the **artificial** \`maxiter=5\` cutoff, reaching maximum scaled 599 residual \`0.0015202525780673786\`. The second call read the exact \`599\`-vector accepted-iterate checkpoint, independently recomputed the same initial residual, rebuilt the finite-difference Jacobian, and accepted three additional full-system steps. It returned a market-root file with all 599 scaled residuals below \`1e-6\`.

| Additional accepted step | Impact \`x_1\` | Full 599 max scaled residual |
|---:|---:|---:|
| 5 — resume input | 0.870026837659 | 0.001520252578067 |
| 6 | 0.868631954353 | 0.0000235839356782 |
| 7 | 0.868649366401 | 0.00000664427261527 |
| 8 — returned root | **0.868649787341** | **0.000000323557477956** |

The returned root is a **nominal depreciation of +16.7595781533 log points**, relative to the same Julia pre-shock terms \`x_pre=1.0271423614420732\`. The source MATLAB guess \`x_1=0.909858511411\` was just an initial guess, not a fixed-price constraint.

The original Stage 5 strictly accounted upper-price saved market root was independently obtained earlier on the owner PC, using the same Julia Euler-polished boom and credit-reversal context:

| Diagnostic | Stage 5 appreciation-side market root | Stage 8 depreciation-side market root |
|---|---:|---:|
| Impact \`x_1\` | **1.156101137755** | **0.868649787341** |
| Nominal log depreciation points | **−11.8272715755** | **+16.7595781533** |
| Max scaled 599 market residual | 1.23028508779e−7 | 3.23557477956e−7 |
| Max asset market residual (level) | 3.64089709493e−6 | 9.57533740120e−6 |
| Max labor market residual (level) | 9.62533446280e−8 | 2.85459880400e−7 |
| Max actual home-goods residual (level) | 6.67465300808e−8 | 1.67805606965e−7 |
| Max positive-mass interior household budget error | 4.44089209850e−16 | 4.44089209850e−16 |
| Max Euler/KKT violation (relative) | 9.99758604272e−6 | 9.98091485724e−6 |
| Positive-mass Euler violation above 1e−4 | 0 | 0 |
| Distribution mass discrepancy | 1.22124532709e−14 | 1.22124532709e−14 |
| Terminal capital gap | 0.010456242789 | **0.004005510218** |
| Terminal interest gap | 4.64494953865e−6 | **2.87337245442e−5** |
| Terminal wage gap | 0.000214417451 | **0.0000821510493** |
| Terminal profits gap | 3.55584496453e−5 | **0.0000136237229** |
| Terminal maximum household budget gap | 0.011342184578 | **0.076587682526** |
| Terminal mass-weighted absolute budget gap | 0.000180930938 | **0.000824716623** |

All values in this table are from actual **user-supplied Julia execution logs**. The files themselves are on the owner's PC and **were not uploaded to GitHub**. No GitHub CI has yet re-evaluated both files against one reconstructed \`Stage3Context\`.

### Economic interpretation: a genuine multiplicity candidate, not yet a theorem

These results constitute strong **numerical evidence for two distinct finite-horizon market-clearing roots of the same corrected 599-equation model**, provided both files pass independent reconstruction against the identical initial state, household operator, external-credit path and terminal continuation. The latter independent file-to-file verification is not yet complete.

Crucially, the \`x_1=0.86865\` solution does **not** require restoring the original MATLAB's coarse-grid consumption-floor/PCHIP budget bug or accepting its original gross-capital accounting violation. Stage 8 computes the lower branch using the corrected strict portfolio and household policy inside every GE residual evaluation.

The two roots yield opposite signed exchange-rate responses to the **same foreign-credit reversal**. If terminal and horizon checks pass, the key economic research problem shifts from reproducing depreciation to **what selects one branch**. Different nonlinear solver seeds alone are not an economic selection mechanism.

## 2. Why full economic certification still fails

Both paths meet the stated \`1e-6\` **finite-period 599-market** tolerance and have small interior budget/Euler discrepancies. However, terminal **household** budget consistency fails the project's independent \`1e-6\` criterion, particularly on the lower-price branch.

In the existing source \`stage3_household_backward!\`, the date-\(T\) household consumption/savings policies are copied from the **calibrated stationary continuation** rather than solved at the actual date-\(T\) prices. The independent diagnostic instead evaluates that copied terminal policy at the actual \`w.rate[T]\`, \`w.wage[T]\`, \`w.profits[T]\`, and \`w.terms[T]\`. Hence the time-\(T\) household budget may fail even while all interior 599 market equations clear.

The terminal gap must be **diagnosed rather than dismissed as truncation**. In particular, the largest residual may come from the far-tail occupied state, and the stationary copied policy may itself have an interpolation/discretization residual at those exact nodes. The decomposition should distinguish those effects before increasing the transition horizon.

For a fixed terminal wealth node \`a\`, productivity \`s\`, and the copied terminal \`(c_ss(a,s),a'_ss(a,s))\`, write the actual date-\(T\) household budget as

\`\`\`text
R_T(a,s) = (1+r_T)*a + (w_T*L_T)*s + profits_T
           - m(x_T)*c_ss(a,s) - a'_ss(a,s).
\`\`\`

Adding and subtracting the calibrated terminal stationary prices gives an **exact identity**:

\`\`\`text
R_T(a,s)
  = R_ss_policy(a,s)
    + (r_T-r_ss)*a
    + (w_T*L_T-w_ss*L_ss)*s
    + (profits_T-profits_ss)
    - (m(x_T)-m(x_ss))*c_ss(a,s),
\`\`\`

where \`R_ss_policy\` is the stationary budget evaluated at the **same node and copied policy**. This is a purely accounting identity. It does not establish causality or permit changing the economic model to make the diagnostic disappear.

The root-comparison audit computes this decomposition at every **positive-mass terminal household state**, including the maximum-error state and the signed probability-mass-weighted contributions. This identifies whether the issue is terminal price/return mismatch amplified by wealth, or an error already present in the stationary policy/discretization. Do not cut or renormalize tail mass.

## 3. Independent same-model comparison (NO new nonlinear GE solve)

Added script:

\`scripts/stage8_compare_saved_roots.jl\`

It reads the owner's two already computed untracked CSVs:

- \`results/economic_reproduction/economic_flexible_leveraged_euler_polished_STRICT_PORTFOLIO_MARKET_ROOT_UNCERTIFIED.csv\`
- \`results/stage8_depreciation/source_initialized_STRICT_MARKET_ROOT_UNCERTIFIED.csv\`

It **reconstructs one and the same Julia \`Stage3Context\`** from the archived pre-shock Euler-polished boom, and reads both saved root paths using the preexisting \`read_economic_root_guess\`, which checks the calendar horizon, inherited \`K_1\`, and EVERY period's foreign-credit path. For each candidate, it independently reruns the entire Euler-polished household backward solver, strictly scaled gross-portfolio revaluation, the full distribution forward operator, and the full **599 market residual**, using the same \`Stage3Context\`. A saved candidate fails if its recomputed maximum scaled residual exceeds \`1e-6\`.

Only **after** that check does the script build a diagnostic \`TransitionSolveResult\` wrapper (zero new root-solving iterations) and call the existing \`stage3_diagnostics\`. It rejects invalid date-1 household affordability or grid support, compares the resulting upper/lower sign relative to the same \`x_pre\`, and calculates the exact terminal budget decomposition above. No code changes are made to the Julia equilibrium equations.

The script writes:

\`results/stage8_depreciation/two_market_root_comparison.csv\`

It explicitly records \`two_distinct_finite_horizon_market_roots=true\` only after both root paths pass the independent residual checks; \`infinite_horizon_multiple_equilibria_certified=false\` remains in the output.

**Run on the user's Windows PC, from \`FMR-2020\`:**

\`\`\`powershell
Set-Location C:\github\22p8z0omwz9xu
git switch research/fmr-stage8-depreciation-continuation
git pull --ff-only
Set-Location FMR-2020
julia --project=. test/test_stage8_two_roots.jl
julia --project=. scripts/stage8_compare_saved_roots.jl
\`\`\`

The CI runs **only algebraic tests and Julia syntax parse**, not this PC-data-dependent actual two-root comparison (the expensive root files are not checked into the repository). Do not report the comparison as independently executed until its actual CSV/log is available.

## 4. Next economic computation

1. Run **the independent two-root comparison once** at T=200. If either path fails the exact same 599 map or inherited state checks, it is not evidence of multiple equilibria in one economy.
2. Use the terminal-budget decomposition to determine whether the current tail and terminal-state price mismatches explain the much larger low-branch max terminal budget gap \`0.07659\`. Confirm the mass at the worst state and all positive-mass terminal nodes.
3. Only then solve an **extended-horizon economy**, preserving the full initial distribution, exact foreign-credit shock timing, strict portfolio, no mass clipping, and branch-consistent warm starts. Since extending \`T\` changes the terminal stationary closure, re-evaluate **both** roots and compare their impact \`x_1\` and date-1 allocations across horizons. A single longer-horizon low root is not enough to confirm same-model multiplicity across horizons.
4. Any claim of full infinite-horizon equilibria or economic selection must await household-terminal and convergence checks. An algorithm landing on one root is not an endogenous selection principle.

**Current decision:** TWO **provisional** finite-horizon strict-portfolio market roots, opposite currency responses, **independent same-context recheck pending**, full infinite-horizon certification pending. No additional multistart or ad hoc solver tuning is justified.
