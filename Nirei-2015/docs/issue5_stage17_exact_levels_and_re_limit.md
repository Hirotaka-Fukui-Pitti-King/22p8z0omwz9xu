# Stage 17 — exact realized household allocation and the remaining RE obstruction

**Source:** Nirei (2015), §§3.1–3.4, Appendix A, equations (7)–(19), (26)–(37).

## State and timing

The current finite-N firm array `(k_i,t, a_i,t, lambda_i)`, the previous marginal utility `mu_(t-1)`, and period-t information on newly revealed next-period productivity `a_i,t+1` are needed to determine firms' next capital and contemporaneous investment. Aggregate `(K_t, mu_(t-1))` alone is NOT a sufficient Markov state for the exact step-function reaction because two different firm profiles can share the same CES capital and differ in mass near the investment thresholds. When prices are predetermined, current price commitments must also be respected; the paper's optimal pricing rule relates these to the next-capital/productivity profile.

The exact finite-N firm problem uses equation (16) and `k_i,t+1=(1-delta)k_i,t lambda_i^kappa` for `kappa in {-1,0,1}`. The already tested `solve_equilibrium!` computes equation (13) capital fixed points, signed investment equation (14), and Section 3.4 ES1 **only when positive-X roots exist**.

## Source-derived nonlinear ex-post static household allocation

Given predetermined `K>0` and selected contemporaneous `X>0`, equations (7),(8),(11),(12), with the paper's KPR utility, imply
```text
Y = K^alpha L^(1-alpha)
C = Y-X
w = c_L*Y/L = psi*zeta*C*L^(zeta-1)/(1-psi*L^zeta)
mu = C^(-sigma)*(1-psi*L^zeta)^(1-sigma)
R = mu_previous/(beta*mu).
```

Eliminating `Y,C,w` yields a **single increasing equation**
```text
X = K^alpha L^(1-alpha)*
    (1+c_L/zeta-c_L/(psi*zeta*L^zeta)).
```
Under `0<alpha<1`, `zeta>1-alpha`, `psi,c_L>0`, and `0<X<K^alpha psi^(-(1-alpha)/zeta)`, there is exactly one interior root, with `C>0` and `1-psi L^zeta>0`. Source implementation: `src/exact_levels.jl` uses bracketed bisection and checks ex-post resource/production/labor/pricing-kernel residuals.

`simulate_exact_level_period!` reuses the unchanged original three-choice firm fixed-point implementation, uses this exact realized allocation, and carries the **realized** marginal utility to the next period instead of forcing `mu_t=mu_bar*(K_(t+1)/Kbar)^h`. This changes the realized macro continuation and may change later firms' predetermined states. Its expected quantities **still** come from Appendix A local continuum forecasts.

## Why full stochastic RE is not solved

1. The published equation (7), `R_t^-1=beta mu_t/mu_(t-1)`, **defines the stochastic discount factor**; the paper does not impose a separate household risk-free asset Euler equation. Equations (9) and (17), rather than a generic invented bond Euler equation, provide the nonlinear **conditional price and investment value expectations** that must hold.
2. Appendix A `E_(t-1)z_t=(B\\D)z_(t-1)` is a local **continuum** forecast, not a conditional expectations operator over exact finite-N firms' full distribution. This implementation retains it and **does not** assert the exact conditional pricing restrictions (9),(17) or complete rational expectations.
3. Appendix C specifies iid Gaussian log productivity, while Section 3.2's one-notch policy assumes bounded innovations. Existing Stage 10–11 work shows strictly positive probability of profiles with **no positive-X capital fixed point** under the unbounded-support law. On those states Section 3.4's `epsilon=log(X/Xe)` is undefined. Hence there is no globally defined original ES1 state transition over the literal shock support without an added restriction or different model. An exact global stochastic RE operator cannot be computed for the **specified** model as-is.
4. Numerical completion or exact ex-post static residuals do not establish global equilibrium existence, nonlinear rational-expectations consistency or Table 2 replication. The new code is an explicitly labeled source-derived **partial/local** improvement; the previous published-model `scripts/run_benchmark.jl` remains unchanged.

## Local-only tests

```powershell
Set-Location C:\github\22p8z0omwz9xu\Nirei-2015
julia --version
julia --project=. scripts/check_julia_syntax.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); using Nirei2015'
julia --project=. scripts/check_published_table1.jl
julia --project=. scripts/check_equilibrium_oracle.jl
julia --project=. scripts/check_exact_levels.jl
julia --project=. scripts/run_exact_level_path.jl 3 48 20261101
```

Run a 20-period N=110000 local test **only** after the exact static regression and small N pass; this is **not** a complete stochastic RE benchmark. Original `scripts/run_benchmark.jl` is untouched. No new CI workflow was created or triggered; commits are marked `[skip ci]`. Julia runtime, syntax, dependency/API, and numerical outcomes remain **unverified in this assistant environment** pending actual Windows results.
