# Stage 3: actual author MAT fixtures and explicit flexible-price goods clearing

## Source-reference export (no nonlinear root solve)

The ZIP `fmr_author_reference_fixtures.zip` is **not tracked** in this repository. It was a user-supplied/generated bundle. Trying to unpack it without first creating or copying it fails. The supported reproducible source is the original MATLAB `.mat` workspaces.

## PowerShell entry point and first saved-root goods test

On Windows **PowerShell**, first update the working branch and run the
checked-in launcher. It derives its paths from the script location, NOT
MATLAB `pwd` (which may remain in OneDrive/Documents/MATLAB even after
passing `matlab -sd`).

```powershell
Set-Location C:\github\22p8z0omwz9xu
git pull --ff-only
Set-Location FMR-2020
powershell -ExecutionPolicy Bypass -File .\scripts\export_author_reference.ps1
```

The launcher checks all three author workspaces, launches MATLAB using
absolute project paths, verifies all twelve expected nonempty CSV files,
and exits on MATLAB failure. No author GE root-finding is repeated.

```powershell
julia --project=. scripts/evaluate_author_explicit_goods.jl results/author_reference
```

The new Julia script reconstructs the shock-date context directly from
the saved MATLAB distribution and inherited scalars, **without first
root-solving the credit expansion**. It evaluates the original unleveraged
and leveraged `K,L,terms` paths on the explicit 599-dimensional goods
map and prints the first-period goods gap, the maximum goods gap, and
the 399 market-block residual. A strict household infeasibility is
reported as such, not mistaken for a successful GE.

To attempt an actual nonlinear 599-variable root find **only after
reviewing the evaluation**, explicitly pass `--solve`:

```powershell
julia --project=. scripts/evaluate_author_explicit_goods.jl results/author_reference --solve
```

This is a potentially expensive numerical solve, not a claimed result.
Its output is accepted and saved only after the implemented
independent asset, labor, goods, budget and probability checks pass.
The diagnostic evaluation is the recommended immediate first run.

## Failure of the author's nested 399-dimensional flexible-price solve

The original `Transition/solve_transition_bis.m` solves the first-period terms of trade with `fsolve(aggr,0.8)` and does **not** check its scalar return residual or exit flag. At the actual saved author 399-vector, it can return the minimizer of a strictly positive goods-market gap. Therefore raw asset and labor residuals under 1e-8 are not proof of a completely closed GE. The unleveraged source goods residual at date 1 is ~0.024336; with 1/16 leverage it is ~0.002110, even though the source 399 asset/labor market conditions are small. Actual household goods-market gaps can be larger due to wealth revaluation/discrete PCHIP policy budgets.

## Separate economically closed 599-dimensional flexible-price system

The new `src/flexible_full_goods_equilibrium.jl` adds `stage3_flexible_full_residual!` and `solve_stage3_flexible_full`. It is **not** a replacement for `solve_stage3_flexible`. It keeps the original flexible-price markup, assets, labor, inherited foreign debt revaluation, EGM, and forward distributions. It makes `terms[1:200]` explicit as additional unknowns and adds 200 *realized household* home-goods clearing conditions. It uses the existing Rotemberg-price routine at exactly zero nominal adjustment cost, which analytically reduces to the flexible-price markup.

For 200 dates the length-599 candidate is ordered as:
`[K[2:200];L[1:200];terms[1:200]]`.
Example only after constructing a **verified** `Stage3Context`:
```julia
x0 = vcat(initial_transition_guess(ctx.economy), fill(1.0, 200))
sol = solve_stage3_flexible_full(ctx; x0=x0, tolerance=1e-6,
                                  maxiter=70, max_evaluations=5000)
```
A new physical capital / labor / terms solution is **not** certified until the returned residuals pass independent household-budget, mass, goods-market and Euler checks. The original fixed-rate 599 equations are a distinct nominal-rigidity regime.

## Validation boundary

- The new MATLAB exporter and 599 GE source code were committed without local MATLAB/Julia execution (neither runtime is present here).
- Automatic GitHub Actions source/Julia tests are triggered by commits; check the latest run before claiming tests pass.
- Do not mistake an accepted trial point, a source-author nonzero scalar goods gap, or merely matching the adapted original 399 map for a closed equilibrium.
- No new 399 or 599 root is reported by this addition, and Figures 1/4/5 remain unverified until the independent GE checks succeed.


## Owner-PC result: original saved MATLAB roots and new fixed audit (2026-11-03)

The original author MATLAB MAT files were successfully converted into 12 CSV fixtures. The author's two flexible solutions satisfy its 399 external asset/labor conditions but do not satisfy the first-period goods equation. The direct Julia 599-equation evaluation found:

| Actual author saved path | Julia scaled max K/L residual | Julia scaled max realized goods residual |
|---|---:|---:|
| Unleveraged flexible | 1.4041195863789824e-6 | 0.007946781106423637 |
| Leveraged flexible 1/16 | 7.637305550045603e-5 | 0.0014468323869371936 |

The source MATLAB fixed 599-root has raw max residual 5.962021987215849e-10, but Julia's raw asset residual differs by as much as **0.17427609901943697 (residual #9)**. The original and Julia date-1 trade terms exactly coincide, while the inherited fine-grid pre-distribution differs by at most 1.3604591926118757e-5. The fixed discrepancy must not be attributed to the missing nested flexible trade root.

To locate the first fixed source-to-Julia mismatch (prices, initial revaluation, EGM, or mass transport), run from PowerShell in the FMR-2020 directory after pulling the branch:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\export_author_reference.ps1 -WithFixedAudit
julia --project=. scripts\diagnose_author_fixed_map.jl results\author_reference
```

This extracts ORIGINAL saved fixed prices, profits, the t=1 shock-date household distribution, realized household savings, and t=1/t=9 fine consumption/savings matrices. It compares each block at exactly the same saved 599-vector. This is NOT root finding.

For a genuine full 599-dimensional flexible goods-clearing solve initialized at the original leveraged saved path, use:

```powershell
julia --project=. scripts\evaluate_author_explicit_goods.jl results\author_reference --solve-leveraged
```

This solve has not converged or been validated in the current record. Unlike audit mode, solve mode requires asset, labor, goods, household budget, mass and Euler checks before accepting a candidate.

## Stage3 computational decision (2026-11-03)

The owner PC has already converged both explicit 599-variable flexible-price **market residual maps** (unleveraged: max scaled 6.239591057958354e-7 in three iterations; leveraged: max scaled 8.778387995206213e-7 in two). Probability mass and household budget identities hold to about 1e-15. Both were rejected by an *additional strict* max relative Euler/KKT cutoff of 1e-4 (observed maxima 3.6509e-4 and 3.8339e-4); weighted mean Euler errors were only ~3.1e-6 and ~3.0e-6. This is **not a nonlinear market solver failure**. The relevant outstanding test is whether these maximum errors are localized to rare fine-grid states or reflect systematic policy interpolation problems. The diagnostics now print the worst occupied (period,asset,state) and the probability share above 1e-4.

The market-clearing root is now **retained even if Euler accuracy remains pending** in `julia_599_MARKET_ROOT_EULER_PENDING_NAME.csv`, and the 599-vector is saved in `julia_599_MARKET_ROOT_NAME.checkpoint` for reuse. Only passing strict Euler checks writes a `julia_full_599_NAME.csv` file. The market root is not described as a fully certified household GE.

MATLAB's original `backsolve_egm.m` produces fine savings decisions only for dates 1 through T-1. Its final-date fine consumption is externally imposed. Hence `dec_tran_fine_fixed` has **199**, not 200, slices. The actual user MATLAB exporter failed because of our incorrect 200-slice assertion; that has been corrected, along with its synthetic test fixture. When rerunning, use the opt-in fixed audit command from the previous section.

**Execution policy:** PC handles actual 399/599 transition and source MAT-file comparisons. GitHub Actions runs quick source/syntax/unit tests; the heavy Octave parity job is now *manual workflow_dispatch only*, not repeated on every push or pull request. No remote CI market root is treated as an actual research result.

## Economically correct finite-horizon stationary-continuation closure (2026-11-03)

**Correction to all older Stage3 outputs:** the previously solved 599 paths used the source MATLAB terminal shortcut `C_T=C_{T-1}` in the labor/goods equations. They are market roots of that old numerical system, **not yet validated equilibria** under an economically coherent infinite-horizon terminal closure. Do NOT cite those residuals as a fully solved economic transition after this change.

An infinite-horizon model truncated at `T=200` has real economic transition dates `t=1,...,T-1`; date T is the start of the externally supplied stationary continuation household policy. It is mathematically inconsistent to impose both `L_T` and `ToT_T` as free transition variables and force the terminal steady household consumption to clear hypothetical final-date labor/goods markets at prices that have not converged to the steady state.

The corrected explicit 599 residual system has **199 asset markets, 199 labor markets, 199 actual home-goods markets, and 2 terminal anchors**:

```text
R_A(t) = (E_mu[a'(t)] - B_t)/q_t - K_{t+1} = 0, t=1,...,199
R_L(t) = L_implied(C_t,w_t) - L_t = 0, t=1,...,199
R_H(t) = Y_t-I_t-cost_t-exports(ToT_t)-C_H,t = 0, t=1,...,199
L_200 = L_calibrated_SS
ToT_200 = ToT_calibrated_SS
```

Unknowns remain `K[2:200]` (199), `L[1:200]` (200), `ToT[1:200]` (200). The terminal `K_200` proximity, stationary household budget feasibility at actual period-200 prices, and Euler/KKT accuracy are **independent truncation/optimality checks**. Date-200 market residuals are not fraudulently set to zero by copying `C_199`. A numerical market root does not by itself certify a valid infinite-horizon equilibrium.

Re-evaluation must be done on the owner PC using the ORIGINAL MATLAB paths **only as initial numerical guesses**, after updating this branch:

```powershell
Set-Location C:\github\22p8z0omwz9xu
git pull --ff-only
Set-Location FMR-2020
julia --project=. scripts\solve_economic_fixed_ge.jl results\author_reference
julia --project=. scripts\evaluate_author_explicit_goods.jl results\author_reference --solve-leveraged
```

These new equations use distinct checkpoint/output filenames with `TERMINAL_BC`; old market roots must never be silently accepted as new full-GE results. Heavy transition calculations remain PC-only. The new terminal-boundary regime requires observing its actual numerical outputs before claiming convergence.


## Economic equation correction and independent Euler test (2026-11-03)

The source MATLAB transition labor condition is **not the same** as
the model's union wage-setting FOC when varphi differs from one.
Paper equation (5) and the Cobb-Douglas consumption calibration
gamma=theta=1 imply

    L_t^varphi = [chi * (epsilon_w - 1)/epsilon_w] *
                 w_t/(chi_dis * C_H,t).

The old MATLAB/Julia transition code instead set L_t equal to the
right-hand side. With varphi=0.5 this changes the marginal labor
response; it also contradicts the stationary union FOC already used
in steady_state.jl. Both Stage2 credit expansion and Stage3 sudden
stop have now been corrected. Any previously saved paths under the
old labor formula are NUMERICAL INITIAL GUESSES ONLY, not verified
solutions of the economic model.

The interpolation-only EGM policies have an observed maximum relative
Euler/KKT discrepancy of about 4e-4 in low-wealth, lowest-productivity
states. src/accuracy.jl now contains an **independent constrained
one-state Euler root** (household_best_response_saving). The Stage3
post-solve audit recomputes the directly optimal savings of the worst
occupied state conditional on the same continuation policy, without
source MATLAB or a consumption floor. It reports the discrepancy
between the EGM-interpolated asset choice and that one-dimensional
Euler/KKT best response. It deliberately does NOT claim that a zero
best-response discrepancy by itself validates the continuation value
function or the full infinite-horizon economic equilibrium.

## Proper horizon-truncation comparison

The existing finite-horizon 3T-1 closure uses asset markets at dates
1:T-1; labor and realized home-goods markets at dates 1:T-1; and
stationary-boundary labor and terms of trade at date T. The last two
boundary conditions are a **finite-horizon asymptotic approximation**,
not a theorem establishing the infinite-horizon solution. Terminal
capital, prices and household-budget consistency must be independently
checked.

Extending the 200-date exogenous quartic credit-reversal path by simply
recomputing its formula with T=300 would make the contraction decay more
slowly, *changing the economic shock*. We therefore preserve the full
original 200-date credit path and append the limiting foreign asset
supply when the simulation horizon is extended. New tail K/L/terms
values are generated for a numerical starting guess only.

The fixed-FX command now supports the original 200-date problem and a
300-date **same-shock** horizon refinement:

```powershell
Set-Location C:\github\22p8z0omwz9xu
git pull --ff-only
Set-Location FMR-2020
julia --project=. scripts\solve_economic_fixed_ge.jl results\author_reference 200
julia --project=. scripts\solve_economic_fixed_ge.jl results\author_reference 300
```

The second command uses the 200-date saved market path as a numerical
guess, and reports the maximum normalized deviation over the initial
40 quarters in K, L, Y, I, w, r, terms of trade and home consumption.
The original credit shock is unchanged. There is no certification
without independent market, Euler/KKT, terminal, and horizon checks.
These calculations are **not yet executed for the revised economic
labor condition**. Do not report the earlier MATLAB-compatible market
roots as final economic results.

Important remaining initial-state issue: the current
scripts/solve_economic_fixed_ge.jl starts from the archived author
credit-expansion shock-date distribution. A final replication must
replace that inherited distribution with one obtained from the newly
corrected Julia Stage2 economic transition, before claiming fully
consistent Figures 1, 4 and 5.


## PC result: corrected economic 599 credit-expansion market root (2026-11-03)

The user's own Windows Julia 1.12.6 run of the corrected union FOC
and explicit full-goods 599-variable **credit-expansion** problem
reported 9 nonlinear iterations and 611 evaluations, maximum scaled
market residual 5.950086127137766e-7, maximum original-units asset
market residual 1.2533693904970278e-7, labor 5.950086127137766e-7,
home goods 3.7728748791820976e-7, household budget 1.33e-16,
and probability mass 1.25e-14.

**This is NOT a certified infinite-horizon equilibrium.**
Independent failures were: maximum Euler/KKT 3.188352638396009e-4
(period 198, productivity state 1), mass above 1e-4 of 0.7248%,
terminal K gap 0.2808580561, and maximum terminal household budget
error 2.4386471598. The independent constrained one-state Euler
best-response chooses savings 0.1244358901, whereas the interpolated
EGM policy chooses 0.1245779020 (difference 1.420118684e-4). Thus
the highest Euler error is a real policy interpolation/optimality error
conditional on the current continuation approximation, not a
miscomputed audit flag.

The user-PC root CSV was saved at
`results/economic_reproduction/economic_credit_expansion_MARKET_ROOT_UNCERTIFIED.csv`.
The new runner never discards this numerical work: it loads K/L/terms
from the saved JULIA CSV, verifies the exogenous B path and inherited K,
recomputes ALL 599 household and market residuals under the current
Julia implementation, and reuses the expensive 599-Jacobian root
*only if the re-evaluated residual is still below tolerance.*
Unsuccessful Euler and terminal checks no longer prevent calculations
of the subsequent credit-reversal market roots, but **all outputs remain
provisional and no final Figures 1, 4, or 5 are generated from failed
economic checks**.

The paper's source initial credit-expansion terminal SS is calibrated
at the date-200 foreign asset stock B_200, NOT the infinite-horizon AR
limit B_infinity. An earlier proposed AR(1) extension to date 300
without recalibrating that SS is wrong (CI caught this mismatch);
the corrected same-experiment horizon refinement holds exogenous
B_t=B_200 after t=200. The full original 200 dates remain identical.

From `FMR-2020/` on the user's Windows PC:

```powershell
julia --project=. scripts\run_economic_replication.jl
julia --project=. scripts\run_economic_replication.jl --horizon=300 --boom-only
```

The **first** command resumes/rechecks the saved 200-date boom and
attempts the unleveraged, leveraged flexible-FX and leveraged fixed-FX
599-variable roots, retaining each for further numerical validation.
The **second** is a separate, much more expensive 899-variable
same-shock boom horizon refinement; it reports the first 40 quarters'
maximum relative differences in K, L, Y, I, w, r, terms of trade and
home consumption. Run it after the 200-date root has been preserved.
Neither command should be called a solved complete replication unless
the independently reported conditions pass. Avoid rerunning a
599-column Jacobian on an already solved boom with unchanged equations.

Stage 3 audit now separately reports maximum Euler errors over the
first 40 dates and the last 20 transition dates, along with
probability-weighted terminal budget errors and terminal wage/interest
deviations. This distinguishes interior economic accuracy from the
numerical truncation layer without hiding either.


## Actual four-root and horizon-convergence evidence, November 3, 2026

The user executed Julia 1.12.6 on their Windows PC after correcting
the non-linear labor union FOC and the Stage3 full-goods residual system.
**All four 599-market systems converged.** These are market roots,
NOT yet households-and-terminal-certified infinite-horizon GE.

| Experiment | Max scaled market residual | Max Euler/KKT | Terminal K gap |
|---|---:|---:|---:|
| Credit expansion 200 | 5.9500861271e-7 | 3.1883526384e-4 | 0.2808580561 |
| Credit expansion 300 | 6.5743535060e-8 | 2.8449041937e-4 | 0.0419838752 |
| No-leverage flexible reversal | 2.3264688896e-7 | 4.1055298560e-4 | 0.0136179332 |
| Leveraged flexible reversal | 1.0497050917e-7 | 3.5860642911e-4 | 0.0106627800 |
| Leveraged fixed reversal | 3.9370336404e-7 | 5.9554357097e-4 | 0.0012221678 |

At identical **200-date exogenous paths**, the first 40 quarters of the
credit expansion differ between T=200 and T=300 by at most
**6.1709093302e-7** normalized, including investment and consumption.
Thus there is no reason to rerun huge horizon extensions before
repairing Euler interpolation. This does NOT prove the terminal steady
household policy is exactly feasible at T=200 or T=300: the T=300
terminal budget maximum is 0.2168190 and probability-weighted
terminal budget absolute error is 0.0030172, both above tolerance.
The large max terminal budget can be dominated by high-asset cells;
the weighted and first-40 measures should be reported separately.

The dominant outstanding **economic numerical error** is nonbinding
Euler optimality at small asset values. Direct one-state conditional
Euler roots differ from existing inverse-EGM decisions by
0.00012 to 0.00024 savings units, in positive-mass cells.
Do not downplay this as merely a printed diagnostic.

### Independent one-pass Euler/KKT interpolation correction audit

To avoid more ~600-column finite-difference GE solves just to
test an interpolation hypothesis, the existing research runner now
supports a **one-pass household-only audit** at the SAVED 200-date
credit-expansion 599-market root:

```powershell
Set-Location C:\github\22p8z0omwz9xu
git pull --ff-only
Set-Location FMR-2020
julia --project=. scripts\run_economic_replication.jl --boom-only --audit-euler-refinement
```

It reloads/re-evaluates the saved Julia 599-market root, calculates
direct constrained Euler/KKT corrections using the *actual next-date
coarse policy* and precomputed PCHIP slopes, first on coarse policies
and then fine-grid policies, recursively from date 199 backward.
It forward-transports the actual household distribution again and
recomputes 599 aggregate market residuals at the UNCHANGED K/L/terms
candidate. It prints how many policies were corrected, maximum
correction in saving, supported-grid post-Euler violation, original
market residual and after-correction market residual. No silent
truncation of desired saving at the upper grid bound.

**Interpretation:** If the post-correction maximum Euler violation is
small and aggregate market residuals are still acceptably small,
the numerical correction is promising. If market residuals shift,
the GE must be solved AGAIN **with the corrected household policy
operator inside every aggregate residual evaluation**. A post-hoc
household correction is never claimed to be a new equilibrium.
The implementation is initially a controlled audit, not the accepted
full equilibrium solution.


## Corrected household-policy operator inside the full GE residual (2026-11-03)

The owner's independent fixed-price household-only experiment on the
stored 200-quarter credit-expansion 599-market root showed:

- 8,079 coarse household policies and 34,099 fine policies corrected
  by solving the conditional Euler/KKT optimal saving choice directly;
- 21 off-grid fine states were not clipped or relabeled optimal;
- maximum Euler violation on the supported correction grid fell to
  9.996579082528452e-6 (from 3.188352638396009e-4);
- largest change in fine-grid saving was 3.241376478071989e-4;
- maximum full 599-market residual increased from
  5.950086127137766e-7 to 3.241547888277463e-5;
- max componentwise market-residual shift was 3.2680082299951785e-5.

**Thus the old GE root is NOT a root of the improved household-policy
system.** A post-hoc household correction cannot be reported as GE.

The fixed-FX and full flexible-FX market residual functions now accept
`polish_households=true`, applying the independent
backward Euler/KKT corrector BETWEEN the base EGM recursion and the
forward distribution. The same corrected micro policy is used on every
nonlinear iteration, Broyden trial step, finite-difference Jacobian
column, and accepted root. No separate/generalized GE solver is created.

Only a targeted 200-date credit-expansion recalculation is authorized
at this stage, from the previously converged 599 aggregate path. This
uses a DISTINCT output file; the source interpolation-only GE roots
remain available and are never overwritten.

```powershell
Set-Location C:\github\22p8z0omwz9xu
git pull --ff-only
Set-Location FMR-2020
julia --project=. scripts\run_economic_replication.jl --boom-only --solve-polished-boom
```

The command must report current market maximum, first-40 and overall
Euler/KKT accuracy, the final-period household budget/asset gaps, and
save its own explicitly unverified economic market-root CSV. The
terminal budget gap, stationarity consistency and asset-grid convergence
remain separate requirements. Do NOT plot or certify Figures 1/4/5
solely because the corrected 599-market residual converges.


## Actual Euler-polished 200-date GE outcome and corrected 300-date test (2026-11-03)

User-PC Julia 1.12.6 run:
`julia --project=. scripts/run_economic_replication.jl --boom-only --solve-polished-boom`

- Independent corrected 599-market equilibrium root converged after 2 nonlinear
  steps / 604 residual evaluations from the saved baseline-EGM market root.
- New maximum scaled market residual: **2.558332579116751e-8**.
- Asset, labor and goods market residuals (native units):
  **2.1377693215640647e-7**, **2.558332579116751e-8**,
  **1.976769303801973e-8**.
- Maximum Euler/KKT: **9.992665155151346e-6** (vs baseline
  3.188352638396009e-4); occupied Euler states exceeding 1e-4: **0**.
- Weighted mean Euler error: **4.523967749583114e-7** (baseline
  3.1532937519699157e-6).
- Household budget: 4.44e-16; probability mass: 1.22e-14.
- First 40 quarters vs baseline-EGM market root, maximum normalized
  difference: **1.744521140878419e-5**, arising in investment.
- **NOT yet infinite-horizon certified**: terminal capital gap remains
  **0.281054172258294**, terminal budget max **2.4473603018791437**,
  probability-weighted terminal budget **0.0293744196445066**.

The economically relevant first-40-quarter response is insensitive to
the Euler refinement on this grid, even though the old *household
optimality* errors were large enough to correct. But Euler accuracy
and horizon approximation are logically independent.

The Julia runner now accepts a **polished 300-date same-shock comparison**:
```powershell
Set-Location C:\github\22p8z0omwz9xu
git pull --ff-only
Set-Location FMR-2020
julia --project=. scripts\run_economic_replication.jl --boom-only --solve-polished-boom --horizon=300
```

It verifies and uses the existing *polished 200-date root* as the
initial-40-quarter comparison baseline; it may reuse the already
solved 300-date BASE-EGM market root as a **numerical starting guess
only**. It then recomputes Euler-polished household choices and
their distribution in EVERY 899-equation GE evaluation and reports
both first-40 economic-path changes and terminal capital/budget
distortions. The 300-date corrected GE must actually converge on the
user PC before any 300-date outcomes can be claimed.

Existing 200-date results are not discarded or silently overwritten.


### Full four-economy Euler-polished re-estimation (next quantitative output)

The corrected 200-date boom GE root is already saved. Rather than
re-running 200-date baseline shock paths, the following mode solves the
complete **four-equilibrium sequence** with the SAME high-accuracy
household policy operator (`polish_households=true`) inside all
599-dimensional market residuals:

```powershell
Set-Location C:\github\22p8z0omwz9xu
git pull --ff-only
Set-Location FMR-2020
julia --project=. scripts\run_economic_replication.jl --solve-polished-all
```

1. Re-evaluate the archived Euler-polished 200-date credit expansion,
   including all consumption and household-distribution objects.
2. Construct the *new* pre-sudden-stop portfolio distribution and
   capital stock from the **Euler-polished boom** (not the old boom).
3. Independently solve the **nonleveraged flexible-FX**, **leveraged
   flexible-FX** and **leveraged fixed-FX** 599-equation market systems,
   each using Euler/KKT-refined household policy throughout.
4. Reuse the stored base-EGM paths of those same horizons **only as
   numerical starting guesses**; verify their B paths and inherited
   capital against the BASE boom to avoid treating stale 599-GE paths
   as solutions in the new polished model.
5. Save new CSVs under distinct names containing
   `_euler_polished_MARKET_ROOT_UNCERTIFIED.csv`, report
   independent Euler/KKT / market / terminal diagnostics and first-40
   aggregate differences versus the old approximations.

No quantitative FX or welfare conclusion is claimed before this
*new* expensive full-economy run is actually executed by the user.
Even if all four market roots converge, the residual terminal budget
and grid/horizon error tests remain separate certification criteria.
The earlier 200-vs-300 Euler-polished boom horizon command is available
independently to check the shared pre-shock path.


## All four Euler-polished GE roots and 200-to-300 refined horizon: actual PC logs (Nov 3, 2026)

The user-PC solved the first 200-date economic shock experiment
and all three 200-date reversal/FX equilibria with the same
high-accuracy household EGM/Euler/KKT operator. The computed GE root
residuals and worst occupied Euler errors were:

| Corrected economy | Max scaled 3T-1 residual | Max Euler/KKT | Terminal K gap | Max terminal budget |
|---|---:|---:|---:|---:|
| 200-date credit expansion | 2.55833e-8 | 9.99267e-6 | .281054 | 2.44736 |
| No-leverage flexible reversal | 3.24100e-7 | 9.98613e-6 | .013469 | .002359 |
| Leveraged flexible reversal | 2.65095e-7 | 9.99738e-6 | .010511 | .011091 |
| Leveraged fixed reversal | 1.10523e-7 | 9.97868e-6 | .001374 | .076837 |
| 300-date credit expansion | 2.35698e-8 | 9.99267e-6 | .042182 | .225235 |

On all four 200-date roots, there were **zero positive-mass Euler states
with violation over 1e-4**. Under the refined household operator,
the initial 40 quarters' paths differ from their baseline-EGM
equilibria by at most 1.74452e-5 for expansion, 9.62928e-6 for
unleveraged flexible reversal, 4.76369e-5 for leveraged flexible
reversal, and 6.86487e-6 for leveraged fixed reversal (scaled
max componentwise deviation). This is evidence against economically
material perturbations in those aggregate paths from the inverse-EGM
interpolation correction **at the current asset-grid resolution**.

Most importantly, **with precisely the same Euler-polished household
algorithm and foreign credit path**, the 200-vs-300 date boom's
first 40 economic quarters differ by at most **6.193045294e-7**
(normalized, worst variable investment). So the reported early response
is insensitive to the lengthening of the temporary foreign-credit
expansion experiment, although the terminal budget gap remains large
by the strict max criterion. This does NOT establish full infinite-horizon
accuracy for the sudden-stop/fixed-FX paths; their own horizon checks
remain unperformed.

### Substantive macroeconomic effects, not more solver diagnostics

The next research deliverable is to answer the paper's economic
questions from the FOUR SAVED corrected market roots without re-solving
the entire 599-equation system. A lightweight script now constructs
the same-date continuing-boom counterfactual and exact pairwise
fixed-vs-flexible and leveraged-vs-unleveraged comparisons:

```powershell
Set-Location C:\github\22p8z0omwz9xu
git pull --ff-only
Set-Location FMR-2020
julia --project=. scripts\report_polished_economic_results.jl
```


The script uses only the actual Julia `economic_*_euler_polished_MARKET_ROOT_UNCERTIFIED.csv` files already on the PC, not MATLAB or rerun policies. It validates that all three reversal regimes share the same 200-date credit shock and shock-date inherited capital, constructs log Cobb-Douglas composite consumption

    log(C_t) = chi*log(C_H,t) + (1-chi)*log(C_F,t) + constant,

reports the direct FIXED minus FLEXIBLE FX impact and 30-quarter
average/cumulative differences, and the leveraged minus unleveraged
consumption effect under floating FX. Separately reports each
reversal relative to the time-matched continuing-boom counterfactual,
not an arbitrary assumed steady-state baseline. The output is
**PROVISIONAL** CSV and Markdown, not certified welfare comparison
and not a claim of complete reproduction of original figures.

This is the quickest meaningful next research step; it requires
reading 4 short CSVs, not another 599/899-column Jacobian.
