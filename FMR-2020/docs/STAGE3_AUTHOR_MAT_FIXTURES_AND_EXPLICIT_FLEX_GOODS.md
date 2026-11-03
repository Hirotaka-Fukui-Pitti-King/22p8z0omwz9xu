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
