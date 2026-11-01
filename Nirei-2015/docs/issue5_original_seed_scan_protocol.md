# Original Nirei (2015) seed sensitivity — local Windows protocol

**Purpose:** Test whether the original published three-choice (`{-1,0,1}`) model, with the published sigma=1.5 calibration, original ES1 and original `simulate_period!` can complete a finite random path for some seeds. This is a **different test** from the old fixed `20261101` failure. There are no fixes to household/macro continuation or firm investment policies and NO shock rejection.

**Commit base:** `research/issue5-original-sigma-correction` / `9cf5701b739ce50e3abcd87015059d329b121261` (NOT the newer Stage17 partial nonlinear levels implementation).

**Local Windows PowerShell:** from `C:\github\22p8z0omwz9xu\Nirei-2015`:

```powershell
julia --project=. scripts/check_julia_syntax.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); using Nirei2015'
julia --project=. scripts/check_published_table1.jl
julia --project=. scripts/check_equilibrium_oracle.jl

# Small-N, 2-seed smoke: logs failures rather than throwing on no-X equilibrium
julia --project=. scripts/scan_original_seeds.jl 2 48 20261102 2

# Fixed ex-ante consecutive range of 20 NEW seeds, 6 periods each
julia --project=. scripts/scan_original_seeds.jl 6 110000 20261102 20
```

The scanner takes `PERIODS FIRMS FIRST_SEED COUNT`, writes a quoted CSV incrementally after **each seed** and final summary TOML into `results/original_seed_scan_T*_N*_start*_count*.{csv,toml}`. Records **both successes and failures**. The original benchmark seed `20261101` is already documented as failing at t=5 and is not repeated in the prospective scan.

If any seed passes six periods, run that specific seed for 20 periods using **this exact unchanged original path**:

```powershell
julia --project=. scripts/scan_original_seeds.jl 20 110000 SEED 1
```

Only for those passing T20, run 1100 periods, preferably with the same scanner to preserve output files:

```powershell
julia --project=. scripts/scan_original_seeds.jl 1100 110000 SEED 1
```

For a T1100 path, the CSV has post-burn-in (first 100 omitted) standard deviations of log Y,C,X,L,K in percent, and residual / trend indicators. Check whether the path stays near the steady state, whether **actual level goods market clearing** holds, and whether moments match Table 2. The original implementation exponentiates a first-order macro model, so exact level clearing is **not guaranteed** even if the path completes.

**Interpretation:** A successful seed demonstrates **finite-path numeric viability conditional on its particular shocks**, NOT an almost-surely well-defined ES1 process over Gaussian support, full stochastic rational expectations, or unbiased replication of published moments. Denominators of every predeclared seed batch and each failure must be retained; never silently discard failed seeds or optimize the seed to fit Table 2.

**Execution status at commit:** Julia is not available in the assistant's environment. The newly added scanner has NOT passed local Julia syntax/import/smoke checks; the user must run these. No GitHub Actions workflows are used.
