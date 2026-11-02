# Stage 3: actual author MAT fixtures and explicit flexible-price goods clearing

## Source-reference export (no nonlinear root solve)

The ZIP `fmr_author_reference_fixtures.zip` is **not tracked** in this repository. It was a user-supplied/generated bundle. Trying to unpack it without first creating or copying it fails. The supported reproducible source is the original MATLAB `.mat` workspaces.

On **Windows PowerShell** (not the MATLAB Command Window), run the following from `FMR-2020/` after confirming `C:\github\HANKSOME\exp_flex.mat`, `exp_flex_lev16.mat`, and `exp_both_41_lev16.mat` exist. This invokes MATLAB as a process; the syntax `addpath(...)` is MATLAB syntax and cannot be entered directly into PowerShell:

```powershell
Set-Location C:\github\22p8z0omwz9xu\FMR-2020
matlab -sd "C:\github\22p8z0omwz9xu\FMR-2020" -batch "addpath(fullfile(pwd,'matlab')); export_author_reference_fixtures('C:\github\HANKSOME',fullfile(pwd,'results','author_reference'))"
```

If `matlab` is not on PATH, replace `matlab` with the installed MATLAB executable path (e.g. `& 'C:\Program Files\MATLAB\<version>\bin\matlab.exe' -batch "..."`), or enter these commands into the **MATLAB Command Window**:

```matlab
cd('C:\github\22p8z0omwz9xu\FMR-2020')
addpath(fullfile(pwd,'matlab'))
export_author_reference_fixtures('C:\github\HANKSOME',fullfile(pwd,'results','author_reference'))
```

The resulting directory must contain **12 CSV files**: 3 inherited/boom files and 3 files per each of 3 saved transition solutions. Check before Julia evaluation:

```powershell
(Get-ChildItem .\results\author_reference\*.csv).Count
Test-Path .\results\author_reference\original_boom_candidate.csv
julia --project=. scripts/evaluate_author_matlab_solution.jl results/author_reference
```

This Julia evaluator ALSO computes its own pre-shock boom GE. It can fail independently if that solve does not converge, even if the CSV bundle is complete. No positive result can be reported without its actual output. **Do not replace source reference data with arbitrary synthetic CSVs.**

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
