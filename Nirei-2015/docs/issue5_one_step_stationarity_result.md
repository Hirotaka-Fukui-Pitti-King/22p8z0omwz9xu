# One-step invariant-distribution test — source-derived analytical result

## Source and scope

Nirei (2015), Section 2.6 Eq. (6) and Proposition 4 (published pp. 964–965), Section 3.3 Eqs. (13)–(16) (pp. 969–970), Appendix C (p. 983). Reference implementation read from the user's `dataverse_files.zip`: Nirei–Scheinkman `setparams.m`, `get_equilibrium.m`, `simulations.m`, `simulations_r.m`.

This is a single, decisive test of the *stationary-gap assumption*, not a proof that the exact finite-N equilibrium paths are stationary. The original finite-N CES root solver and ES1 are unchanged.

## Exact gap transition under the original three-choice policy

Define `b_i,t = b_i(a_i,t)`, `s_i,t=log(k_i,t/(b_i,t K_t))/log(lambda_i)`. Original firm law implies

```text
s_i,t+1 = s_i,t
        + [log(1-delta)+log(b_i,t/b_i,t+1)+log(K_t/K_t+1)]/log(lambda_i)
        + kappa_i,t+1,     kappa ∈ {-1,0,+1}.
```

For constant productivity and K, let `d_i=-log(1-delta)/log(lambda_i)`. If `lambda_i(1-delta)>1` then `0<d_i<1`, and the original one-notch threshold rule gives `s_i,t+1=(s_i,t-d_i) mod 1`. Because the rotation preserves Lebesgue measure on [0,1), the conditional uniform gap distribution is invariant in this deterministic environment. No code correction is warranted there.

For IID Appendix C log productivity, `b(a)∝a^q` with **`q=rho/[alpha(1-rho)]=10`** under published `alpha=.26, eta=11, delta=.02`. Holding aggregate K fixed at Kbar for this *partial-equilibrium one-step experiment*, put

```text
r = s-d + Z,         s ~ U[0,1),
Z = q*(log(a_old)-log(a_new))/log(lambda)
  ~ N(0, [q*sqrt(2)*sigma_log_a/log(lambda)]^2).
```

The **original** restricted policy chooses `kappa=+1` for `r<0`, `kappa=0` for `0<=r<1`, and `kappa=-1` for `r>=1`, so `s_next=r+kappa`. Hence `s_next<0` if `r<-1`, and `s_next>=1` if `r>=2`. These have positive probability under an unbounded Normal, even if the original gap was exactly uniform. Therefore the original *uniform distribution supported on [0,1)* is **not an exact invariant distribution** under the published stochastic calibration and restricted one-notch choice set. This conclusion requires no simulation or Julia execution.

### Exact probability, integrating out the uniform gap and λ

Conditional on λ, define `v=q*sqrt(2)*sigma_log_a/log(lambda)`, `d=-log(1-delta)/log(lambda)`, `H(u;v)=u*Phi(u/v)+v*phi(u/v)`. Then

```text
Pr(s_next<0 | lambda) = H(d-1;v) - H(d-2;v).
```

The corresponding upper exit probability follows from the Gaussian upper tail above `1+d`; it is negligible but strictly positive. Quadrature under `lambda ~ Normal(1.028,.004) truncated to [1.02,1.036]` gives:

| Metric | Analytically integrated probability |
| --- | ---: |
| λ(1−δ)≤1 (violation of Section 2.6 support) | **0.639140832%** |
| 1-period s_next<0 with IID productivity | **2.968281375%** |
| 1-period s_next≥1 | **~0.000000001%** |
| Total outside [0,1) | **2.968281377%** |

For λ=1.028: `d=0.7315801`, standard deviation of productivity threshold displacement `v=0.2560574`, `Pr(exit | lambda)=1.94432%`.

This is **a contradiction between Section 2.6's support assumption and Appendix C's calibrated unbounded shocks**, not by itself proof of a transcription bug in `src/firms.jl`. The code's distributions reproduce Appendix C. Changing λ or truncating productivity, or permitting more investment notches, would modify the specified published model.

## Independent numerical experiment (Python, not Julia)

Computed independently from the published policy, 2,000,000 IID sampled firms (NumPy RNG seed 20261101, SciPy truncated Gaussian λ; this does not use Julia's Xoshiro streams). Fix threshold aggregate K to the theoretical Kbar (normalized to 1). Initialize `k0=b(a0)*Kbar*lambda^s`. Compare **identical** initial firm draws with a no-productivity-change control and an IID productivity change; compute candidate next-period CES K and signed X. **This is NOT a self-consistent aggregate-capital fixed point or ES1 equilibrium.**

| Measure | No productivity shock | Independent Appendix C shock |
| --- | ---: | ---: |
| `K1/K0-1` | **−0.000534%** | **−0.073518%** |
| Candidate `X` in Kbar-normalized units | 0.01999489 | 0.01926415 |
| Outside [0,1) | 0.0057% (due to λ support violations) | **2.96895%** |

Baseline Monte Carlo `K0/Kbar=1.0000011088`. The within-draw capital difference due to the productivity shock is `K1_shock-K1_control=-0.00072984` in Kbar-normalized units. This **isolates the one-period failure of the stationary approximation** and suggests a channel for negative investment drift, but does **not** establish that it causes the observed t=5 full-GE collapse, since endogenous K/ES1 is not solved in this experiment.

### Reproduce on the user's actual original Julia code

The committed `scripts/check_stationary_one_step.jl` uses the existing Julia `initialize_firms`, `threshold_coefficient`, `draw_productivity!`, and `aggregate_capital`. It holds the threshold's target K at Kbar, compares matched no-shock and shock outcomes, and verifies the exact gap transition algebra. It does not touch model source, stochastic draws, 3-choice policy, ES1, or the fixed-point enumerator.

```powershell
Set-Location C:\github\22p8z0omwz9xu
git status --short
git fetch origin
git switch research/issue5-stage18-original-implementation-repair
git pull --ff-only
Set-Location .\Nirei-2015
julia --project=. scripts/check_julia_syntax.jl
julia --project=. -e 'using Nirei2015'
julia --project=. scripts/check_stationary_one_step.jl 110000 20261101
```

**Verification:** Standalone mathematical derivation and independent Python one-period numerical computation completed. Julia is **not installed in the assistant's runtime**; newly committed Julia source is **syntax unverified, package execution unverified, and smoke-test unverified** until the user runs the commands locally.

## Implementation decision

- **Initial `initialize_firms`:** Samples `s0~U[0,1)` as Appendix C explicitly says, and normalizes E[a^q]=1. Not a demonstrated Julia implementation error. Finite-N realized `K0!=Kbar` is ordinary sampling unless an exact alternative finite-N normalization is specified; rescaling real physical K would change actual gaps and decisions.
- **Firm transition:** Original code implements precisely `k_next=(1-delta)k_old*lambda^kappa` with the original one-notch thresholds. Not a demonstrated typo or off-by-one bug.
- **Invariant distribution:** **FAIL under the original Appendix C stochastic calibration in a fixed aggregate-K one-period benchmark**, by strict analytical support argument and explicit quantitative test.
- **Source modification required?:** No uniquely source-justified correction while preserving simultaneously (i) exact published 3-choice policy, (ii) Appendix C Normal shocks, and (iii) exact uniform stationarity. These requirements are mathematically inconsistent. No unlicensed production fix is made.
- **Dynamic RE and original T6/T20/T1100/Table 2:** No new run. Previous seed20261101 failed at t5, and 20 predeclared seeds failed before T20. Cannot claim full-model recovery.

The successful result is a **definitive one-period failure of the stationary-distribution assumption for the published calibration**, not a successful replication and not a proof of the sole cause of capital collapse. Issue #5 remains open.
