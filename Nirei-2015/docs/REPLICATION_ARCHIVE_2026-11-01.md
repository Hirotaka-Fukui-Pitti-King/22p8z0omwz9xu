# Nirei (2015) replication — final archival handoff (2026-11-01)

**Status: PAUSED / NOT REPLICATED.** The user requested a halt to further numerical debugging. No 2015 Table 2 replication or global finite-N stochastic rational-expectations equilibrium has been established. This is a record of *negative and partial results*, not a validated reproduction.

## How this archive was made

The initial `main` base was `26868a8ec6786ff0a872ce24764d07dd32ffa2a4`. The authoritative final working tree was `research/issue5-stage18-original-implementation-repair` at `a5f159ca1b16545a1b2a31e3d3540bc9fc56ae15`, including the Stage 19 two-choice first-order continuation and Stage 20 CES-accounting diagnostic. The final archival merge commit has the original main and **all eleven work-branch tips as parents**, keeping all investigations in the main commit DAG. Unique source/data from divergent branches were copied into the archival tree; overlapping experimental changes were *not* silently substituted for the final Stage 18 code. The original 2015 three-action implementation remains intact.

Archived branch heads (Git SHA remains permanently reachable through the main merge ancestry):

| Branch | Final SHA | Distinct investigation |
|---|---|---|
| `research/issue5-author-epsilon-protocol` | `8a7f9d63ed7f42cca65a0430759c5d431c2c27ac` | author-style epsilon protocols |
| `research/issue5-historical-binary-protocol` | `d8bfcff155dea362080d8395ec6bb64e05b78f14` | historical two-choice protocol |
| `research/issue5-original-seed-sensitivity` | `f8a6b879bdf16c168ff429eac863300fa4f35e7d` | original seed sensitivity |
| `research/issue5-original-sigma-correction` | `9cf5701b739ce50e3abcd87015059d329b121261` | published sigma calibration correction |
| `research/issue5-stage8-investment-wedge` | `ee444c13a525ba3699480698471cb31b0d1b7318` | investment wedge |
| `research/issue5-stage10-fixed-point-feasibility` | `3ee1956ac9b436cbaa15c84a74d1c8d7ff1166f8` | finite equilibrium feasibility |
| `research/issue5-stage11-global-equilibrium-domain` | `cc870db10187a3d02378c645a9a15f4573812353` | global equilibrium domain |
| `research/issue5-stage12-multistep-equilibrium` | `d351c52342c95704cd79afea3d8dbb0208a352c2` | multistep equilibrium enumeration |
| `research/issue5-stage13-dynamic-multistep` | `eb76c3449e6aa308ad9829dc2a37b9a6e0e77510` | dynamic multistep |
| `research/issue5-stage17-dynamic-re-equilibrium` | `6174d1f83f45fd556df09321c00f3e8bb3c425c7` | nonlinear exact-level diagnostic |
| `research/issue5-stage18-original-implementation-repair` | `a5f159ca1b16545a1b2a31e3d3540bc9fc56ae15` | latest three-choice/two-choice and Stage 20 |

## Definitive observed outcomes

1. **Canonical three-action** `κ∈{-1,0,1}`: exact finite-N fixed-point enumeration with ES1 and 110,000 firms; original seed `20261101` reaches a state at `t=5` with no admissible strictly positive aggregate investment equilibrium under the implemented log-investment ES1 domain. This does **not** establish absence of all finite-firm capital fixed points. No 1,100-period benchmark or ten-replication Table 2 result was produced.
2. **Book-motivated two-action** `κ∈{0,1}`: the small-`N=8` oracle enumerating all `2^8` policies passed. A 110,000-firm/1,100-period Julia 1.12.6 run using first-order stable-manifold macro continuation completed, but `K` fell from approximately 8.65 to `0.0028847244615743776` and `X_{1099}=3.636627996982054e-5`. This is **not** the original 2015 three-action benchmark.
3. Two-choice post-burn-in 1,000-period moments: SD(log Y) 133.7918%, SD(log C) 119.5518%, SD(log X) 210.2470%, SD(log L) 106.9862%, SD(log K) 210.0951%. Correlations Corr(Y,C) 0.999992, Corr(Y,X) 0.999922, Corr(Y,L) 0.999988, Corr(Y,K) 0.999976. Nirei (2015) Table 2 values respectively 2.28%, 0.86%, 6.60%, 3.06%, 0.26%, 0.822, 0.979, 0.795, -0.018.
4. The previous hybrid marginal-utility continuation reached `t=139` in the two-choice economy and failed at `t=140` with `K_prior=3.1387408929990417` and `mu_prior=4.868312801534069e11` in the *ex-post nonlinear macro residual check*. This is not a proof of nonexistence of a capital fixed point. The first-order stable-manifold closure removed that numerical blow-up but did not restore stationarity. It is a **local approximation, not a complete global rational-expectations solution**.
5. Stage 20 verified exact CES growth decomposition, 110,000 firms × 20 periods, seed `20261101`. Average log-capital growth = **-0.007165915824613479 per period**, decomposed into depreciation `-0.020202707317519445`, productivity reweighting `-0.0000067060113723357935`, and discrete capital adjustment `+0.013043497504278306`. The accounting identity and the proportional scaling of all two-choice capital roots were independently tested. **This explains the arithmetic of capital decline, not the economic reason why the equilibrium selects insufficient investment.**
6. Stage 20 first period: `Knext=8.590419357387953`, `X=0.11305740191076478`, 17 positive-X roots, 53,576 investing firms. The expected-investment forecast was approximately 0.173; the finite-firm conditional investment response at expected K is lower, and ES1 selects a still lower realized investment. Do not equate this wedge mechanically with an ES1-only effect.
7. The capital-gap distribution is near uniform in moments but not strictly within [0,1) for all firms. This is **not sufficient** to establish a stationary aggregate-capital process.
8. **Benchmark calibration:** the actual `src/calibration.jl` default and published Table 1 test use `sigma=1.5`. Older `ALGORITHM.md`/`policy.json` text says 15; treat that as outdated documentation, NOT a specification to change the numerical benchmark. See `scripts/check_published_table1.jl`.

## Verified commands and evidence

- Stage 19 20-period short CI **SUCCESS**: https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/actions/runs/41304261106
- Stage 20 20-period CES accounting CI **SUCCESS**: https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/actions/runs/41305594102
- User's Windows Julia 1.12.6 completed: `julia --project=. scripts/run_two_choice.jl 1100 110000 20261101 1 linear`; local time-series/summary `results/two_choice_rep01_linear.csv` and `results/two_choice_rep01_linear_summary.toml` were created **on the user's Windows workstation**. They are not guaranteed to have been committed to GitHub. Never invent a committed file path for these.
- Stage 20 can be reproduced using `julia --project=. scripts/stage20_capital_accounting.jl 20 110000 20261101` or `--smoke`. CI artifact: https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/actions/runs/41305594102
- Earlier original-three-choice and experimental CI runs **failed or were inconclusive**; see historical GitHub Actions and Issue #5. Do not reclassify them as successes.

## Source and method preservation

- Canonical Nirei (2015): `src/{calibration,expectations,firms,equilibrium,simulation}.jl`, plus `scripts/run_benchmark.jl` and `docs/issue5_stage18_source_to_code_audit.md`.
- Historical investigations: `docs/issue5_*.md`, `scripts/diagnostic_*.jl`, `results/diagnostic_*.csv/.toml`; seed sensitivity files and epsilon protocol scripts have been brought forward.
- The Stage 17 exact-level auxiliary code is archived, not used to overwrite or silently amend the Stage 18 package. Historical source and its module variant remain at the archived Stage 17 commit.
- Historical alternative CI workflows are archived outside `.github/workflows/` so they do **not** schedule superseded experiments when this repository changes. Their original paths remain reachable in their ancestor commits.
- Current active CI may report failure for the unresolved original three-choice simulation. The archival merge must not be represented as proof that the benchmark is passing.

## Important cautions for future resumption

**Do not restart the same 3-choice/2-choice ES1 scans or the 1,100-period two-choice run without a changed, source-supported hypothesis.** Prior work already established (a) the exact finite-N fixed-point oracle for the implemented choices, (b) original t=5 positive-X failure, (c) the 2-choice macro continuation fix but enormous stationary-moment mismatch, and (d) the CES arithmetic decomposition.

The unresolved question is why the original 2015 stationary continuum-based investment forecast differs from achievable finite-firm equilibrium investment at the initial state, including the mapping from exact CES capital, firm policy, information sets, and the approximation used by the author. Root cause and defensible repair **remain unproven**. Any proposed repair needs an independent economically justified implication or original author implementation evidence, not a new arbitrary normalization, seed rejection, shock adjustment, or equilibrium-selection change.

This archive is intentionally a **paused negative-result handoff**. Issue #5 and #6 should be closed as **not planned / paused**, not completed. Reference this file and the actual recorded run logs when reopening a future attempt.
