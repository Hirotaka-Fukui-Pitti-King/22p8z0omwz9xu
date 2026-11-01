# Explicit three-choice negative-investment equilibrium completion (Issue #5)

## Economic reason and exact scope

**Nirei (2015), §2.3 Definition 1, p.959:** choose among **all** finite-N capital fixed points the one minimizing `abs(log(K_next/Kstar))`, which is well-defined whenever capital is strictly positive, including when aggregate **net investment is negative**.

**Nirei (2015), §3.4 p.971:** chooses the solution minimizing `abs(log(X_realized/X_expected))`. Unlike Definition 1, this logarithmic objective is **undefined for X<=0**. The paper offers no explicit continuation rule when every capital fixed point has nonpositive X. The original Julia solver had incorrectly described the resulting failure as lack of a feasible equilibrium; in fact the three capital fixed points remain and satisfy the original three-choice firm policy and exact CES.

### Implemented completion — explicitly opt in

`selection_mode=:original` remains the default and preserves the original positive-X Section 3.4 selector. With `selection_mode=:capital_fallback`:

1. Enumerate **all** 3-choice CES fixed points unchanged.
2. If any have X>0, apply original Section 3.4 `min abs(log(X/X_expected))`, unchanged.
3. Only if none has X>0, choose the all-root closest-capital solution `min abs(log(K/Kstar))` using Section 2.3 Definition 1, where `Kstar=(1-delta)Kprior+(A^(rho-1)/rho)*X_expected`.
4. Keep the **actual signed X**. Reaggregate firm capital/investment independently, check all original numerical policy residuals, and advance (k,a) only after successful firm-level checks. **No** clipping or investment floor, change to λ, shock rejection, ES2 or multi-notch adjustment.
5. Since §3.4's log-investment shock and `log(X)` do not exist for X<=0, to continue the **first-order approximation only**, the resulting `EquilibriumResult.epsilon` stores `log(K/Kstar)/delta` as a *capital-innovation proxy*, NOT §3.4's log-X shock. Contemporaneous Appendix A linearized macro equations use `(X-Xbar)/Xbar` as a *first-order Taylor substitute for log(X/Xbar)* whenever X<=0. The PeriodRecord's literal `logX` is **NaN** for negative X. All these fields are *not an exact nonlinear GE or RE solution* and cannot reproduce Table 2 without an independent signed-X compatible structural closure.

**This hybrid fallback is a new model-completion choice, NOT literally published in §3.4.** It preserves the original model on its original positive-X domain. The distinct §2.3 ES1 capital formulation motivates the fallback but does not mathematically identify it uniquely.

### Prior t=5 expected firm-level result

Previously verified original seed=20261101, N=110000 fixed points at t5:

| K | signed X |
| ---: | ---: |
| 7.922286680029842 | -0.05758656336226558 |
| 7.922473044264564 | -0.0574001973735959 |
| **7.922475106160655** | **-0.057398135458978324** |

Forecast Kstar ≈ 8.146420155827748, so `capital_fallback` should choose the **third** point with largest K, conditional on exactly the same RNG state. The original `selection_mode=:original` must still stop because its §3.4 objective is undefined. The fallback removes the **selection-domain failure** but is not evidence the subsequent macro solution remains economically admissible.

### Code

- `src/equilibrium.jl`: explicit `selection_mode` keyword for `select_es1` and `solve_equilibrium!`; strict default, capital fallback only when no positive-X candidate; signed physical investment preserved.
- `src/simulation.jl`: mode propagation and clearly marked first-order signed-X macro extrapolation; literal `logX` recorded as NaN when X<=0.
- `scripts/run_benchmark.jl`: optional 4th argument `capital_fallback`, with distinct output files `results/simulation_rep01_capital_fallback.csv` and `results/simulation_summary_capital_fallback.toml`; status counts nonpositive selected investments. Original mode and its file paths untouched.
- `scripts/check_capital_fallback.jl`: small, deterministic selection regression. It checks strict all-negative rejection, capital-based choice of nonpositive root, and bitwise same original selection when positive-X roots exist.

### Required Windows Julia 1.12.6 verification (not executed in assistant container)

```powershell
Set-Location C:\github\22p8z0omwz9xu
git status --short
git fetch origin
git switch research/issue5-stage18-original-implementation-repair
git pull --ff-only
Set-Location .\Nirei-2015
julia --project=. scripts/check_julia_syntax.jl
julia --project=. -e 'using Nirei2015'
julia --project=. scripts/check_published_table1.jl
julia --project=. scripts/check_equilibrium_oracle.jl
julia --project=. scripts/check_capital_fallback.jl
julia --project=. scripts/run_benchmark.jl 6 110000 20261101 capital_fallback
```

**Gates:** Only if the six-period signed-X variant completes, attempt `20` periods. Do not run 1100 until 20 passes, and do not label any statistics with negative investment as Table 2 log-investment moments.

**Verification status:** GitHub file writes were made but **Julia syntax parsing, package load, regression execution and benchmark are not verified** in the assistant container (no Julia installed). The user must execute local tests and provide logs. No GitHub Actions should run heavy computation. Issue #5 remains OPEN.

## Model choices not made

No changes to the **published three investment choices** `κ∈{-1,0,+1}`; no change to productivity/λ distributions, initial firms, aggregate capital fixed-point solver, investment floors, or seeds. The new capital rule is an explicit **opt-in fallback**, not a silent rewrite of §3.4.
