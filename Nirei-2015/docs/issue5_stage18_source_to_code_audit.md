# Stage 18 source-to-code audit — original Nirei (2015) implementation

**Source:** Nirei (2015), *Theoretical Economics* 10(3), 953–985; Sections 3.1–3.5, Appendix A (Proposition 5), Appendix B (equations 9,18), Appendix C. DOI 10.3982/TE1611.

**Base code:** `research/issue5-original-sigma-correction` at commit `9cf5701b739ce50e3abcd87015059d329b121261`. This audit deliberately does NOT import Stage17's experimental closure or the seed-screening branch. **No production source changed.** Do not claim any numerical benchmark has been rerun.

## Equation-by-equation verdict

| Paper | Source / function | Observation | Verdict |
|---|---|---|---|
| (10), (16), threshold `k_i,t+1^* = b_i,t+1 K_t+1` | `firms.jl:threshold_coefficient`, `equilibrium.jl:enumerate_fixed_points!` | Productivity `a_i,t+1`; `b=A[a^(rho/alpha)(lambda^rho-1)/(lambda-1)]^(1/(1-rho))`, switches at `d/(lambda*b)` and `d/b`. | Correct relative to published equations. Independent all-roots oracle passes; does not itself verify intertemporal equilibrium. |
| (13) and (15), CES capital | `firms.jl:aggregate_capital`, `equilibrium.jl` | Code uses `a_i,t+1`, matching (15), price-setting, and the explicit definition of `K_t+1`, while printed (13) shows `a_i,t`. | Paper-internal time-index discrepancy, not proven code bug. |
| (14), signed investment | `equilibrium.jl:solve_equilibrium!` | Direct summation of `[(lambda^kappa-1)(1-delta)k_i,t]/N`. | Correct. |
| Section 3.4, ES1 | `equilibrium.jl:select_es1` | Selects min absolute log investment surprise only among positive-X roots; stops when none. | Correct as a *partial* numerical rule. The paper does not specify how to extend log-ES1 to X<=0. |
| (26)–(37), conditional log-linear expectations | `expectations.jl` | Computes `F=B\\D` and checks matrix and linearized restrictions. | Algebra verified; finite-N rational forecast consistency not established. |
| Realized log-linear (27),(28),(30),(33) | `simulation.jl:_realized_macroeconomy` | Exponentiates a solution of the log-linear equations. The paper explicitly says to *approximate* the expectation system in the first order, but exact (12) is `Y=C+X`. | **Approximation-domain ambiguity / exact-level equilibrium failure.** User observed max ex-post level resources gap 12.6%–49.5%; not floating error. This is NOT by itself proof of an incorrect implementation of the *published first-order numerical approximation*. |
| Realized continuation `logmu=h*log(K_next/Kbar)` | `simulation.jl:simulate_period!` | `h` from stable eigenvector of Appendix A conditional 2-state linear system. Applies it to exact finite-N realized capital. | Additional bounded-path closure; compatible with a *local* stable-path approximation, but nonlinear conditional expectations (9),(17) are not verified. No unique source-derived nonlinear replacement established. |
| Initial gap distribution, Appendix C | `firms.jl:initialize_firms` | iid uniform `s_i,0`, iid lognormal `a_i,0`, truncated-normal permanent `lambda_i`; finite-N capital is computed rather than forced to continuum Kbar. | Follows explicit Appendix C draws. Whether joint distribution is exactly stationary under finite-N endogenous dynamics is NOT established. |
| Section 3.2 lambda support vs Appendix C | `calibration.jl` | Theory requires lambda>(1-delta)^(-1)=1.02040816..., whereas Appendix C's truncated Normal with mean 1.028, SD .004 and 2-SD bound includes [1.020,1.02040816]. | **Internal paper specification conflict**, not code error. Under the stated truncated Normal, approx 0.639% (=703 of 110000 in expectation) of drawn lambdas violate the theory support restriction. Do NOT truncate or rescale without labeling a model change. |
| Section 3.2 bounded productivity vs Appendix C | `firms.jl:draw_productivity!` | Theory assumes bounded shock support; Appendix C specifies normal log productivity (unbounded). | **Internal paper specification conflict**, not code error. Prior Stage 10–11 establishes possible paths with undefined positive-X log-ES1. |

## Why Stage17 is not an exact general equilibrium repair

`research/issue5-stage17-dynamic-re-equilibrium/src/exact_levels.jl` imposes `w=c_L Y/L` on **realized** variables. In the published model, firms set prices **one period in advance** and labor adjusts passively to demand (Section 3.2, Appendix B derivation of Eq. 9); the constant labor share is a **steady-state / conditional optimal-pricing** restriction. It is not a general exact *realized* factor-price first-order condition under predetermined goods prices. Therefore the Stage17 static residual PASS does not certify the published model's complete realized equilibrium.

Even if `mu` is specified exogenously by the linear stable eigenvector, the exact conditions (7), (8), (11), (12), and the KPR marginal utility relation alone do NOT yield a globally unique realized allocation. Eliminating Y,C,w yields:
```text
mu = [K^alpha * L^(1-alpha) - X]^(-sigma)
     * (1 - psi*L^zeta)^(1-sigma)
Y  = K^alpha * L^(1-alpha),   C=Y-X
w  = psi*zeta*C*L^(zeta-1)/(1-psi*L^zeta).
```
For sigma>1, the function mu(L) diverges both near C=0 and near the leisure domain boundary; for some levels there are TWO interior labor solutions. An internally consistent numerical example using the published calibration and the already observed initial firm K/X is:
```text
K=8.6503730193, X=0.0693682701, mu=1.3321324719
L≈0.4417553655, C≈0.8879660173, Y≈0.9573342874, w≈0.9747472988
L≈0.9268657905, C≈1.5872464977, Y≈1.6566147678, w≈20.8794564279
```
Both satisfy these exact static equations with zero level resources discrepancy. They are **not both certified dynamic equilibria**: future expected prices and other dynamic restrictions may exclude either/both. The example demonstrates why adding `w=c_L Y/L` without the correct time-t predetermined price system is *not* a valid proof of an exact nonlinear solution. The scalar counterexample was independently evaluated with algebraic numerical calculations, **not Julia**, and is not a new benchmark run.

## What is, and is not, identified

1. **Confirmed:** Original log-linear reporting has large **realized level accounting residuals** on short surviving paths. These are first-order approximation errors far outside an infinitesimal neighborhood, not residuals of a solved exact level GE.
2. **Confirmed:** The historical original 20-seed T20 test yielded **0/20** completed paths, with failed prefixes ending in no positive-X fixed point. This observation does not imply no successful seed exists anywhere.
3. **Confirmed:** The Stage17 alternative makes realized static levels satisfy its imposed equations, yet fails at t5; its added realized constant labor share is not justified by the source.
4. **Not established:** A *specific incorrect* threshold, index, RNG, CES aggregation, or 2x2 Appendix A matrix term that can be corrected with a paper-backed replacement and make T20/T1100 succeed.
5. **Not established:** A fully specified stochastic conditional equilibrium operator that connects exact finite-N firm distribution to (9),(17), and a globally defined log-ES1 over the original support.

**Decision:** Do not commit a speculative production patch to `src/`, do not retune seeds or productivity supports, do not replace investment with floors or change ES1, and do not rerun T20/T1100 on unchanged code. A valid source-level repair requires a proved **wrong expression or timing** together with the source-derived replacement and an independent before/after test. A first-order approximation is not mechanically a coding error; the numerical violation must be classified accurately.

**Status:** Stage18 documentation-only audit. New code: none. New Julia syntax tests: not executed; benchmark: not executed. Previous Windows source branch Julia 1.12.6 syntax/calibration/fixed-point oracle passed. Issue #5 remains OPEN.
