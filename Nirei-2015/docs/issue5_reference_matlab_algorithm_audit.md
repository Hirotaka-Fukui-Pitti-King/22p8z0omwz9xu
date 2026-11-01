# Actual reference-code audit: Nirei–Scheinkman *Repricing Avalanches* vs Nirei (2015)

Date: 2026-11-01. Working branch: `research/issue5-stage18-original-implementation-repair`, based on the original corrected 2015 three-choice implementation, not Stage17.

## Verified source artifacts

Read **all** five requested files from the user's uploaded `dataverse_files.zip` in full: `README.rtf` (48 lines), `setparams.m` (31 lines), `get_equilibrium.m` (70 lines), `simulations.m` (141 lines), `simulations_r.m` (217 lines). These are the actual author MATLAB files; do not confuse their pricing model with the 2015 investment model.

Also directly read `vgregory757/KhanThomas/Code/KhanThomas.jl` (not official 2008 authors' code) and `MCHatcher/OBC-multiple-equilibria/Select_solution.m`, `OBC_Example_2_stoch_sims.m`, `README.md` from GitHub. Did not inspect source files in the official Winberry replication archive; no algorithm-level claims are attributed to that code.

Reviewed original 2015 published paper pp. 964–966, 968–972, 980–983 and original corrected Julia sources `src/calibration.jl`, `src/firms.jl`, `src/equilibrium.jl`, `src/expectations.jl`, `src/simulation.jl`.

## Reference-code algorithm, with concrete lines

**Stationary:** `setparams.m:1–31` defines 40,000 firms, 2 productivity states, menu-cost and Calvo rates. `get_equilibrium.m:9–23` performs *outer bisection on steady wage w*, calling `fzero` for two within-state repricing threshold ratios q_a; then computes lower/reset/upper price boundaries and stationary aggregates `ss` (`get_equilibrium.m:24–70`). This is a nested equilibrium computation **for the pricing model**, not a required replacement of 2015's closed-form continuous-firm steady state.

**Sampling and exact normalization:** `simulations.m:60–72` samples `n/2` prices conditional on each of two productivity states by the **conditional stationary inverse CDF**, then divides all sampled prices by their exact empirical CES aggregate `mean(p^(1-eta))^(1/(1-eta))` to impose price-index unit normalization. It does **not** sample a uniform conditional price gap. The cross-sectional distribution is explicitly conditional on the productivity state.

**Propagation:** `simulations.m:75–94`: after one Calvo-induced reset, recompute `ePnew=mean(p_i^(1-eta))`; the implied change in aggregate P triggers additional *downward threshold crossings*, and `while flag>0` repeats their adjustment and reaggregates. `simulations.m:95–130` treats upward and downward price changes differently following a productivity switch. `simulations.m:132–136` re-normalizes the relative-price profile to the new aggregate P and records the event. This is a **shock-triggered monotone cascade conditional on the initial event and direction**; it is not enumeration of all mutually consistent equilibria under the 2015 ES1.

**Extended wage feedback:** `simulations_r.m:89–137` adds an additional wage-induced drift, a first passage to the threshold and corresponding price cascade *before* the next Calvo shock, then `simulations_r.m:140–212` repeats the Calvo response and updates the relative-price profile to account for wage changes. Wage adjustment is an additional endogenous feedback in *that model*, not automatically a valid 2015 household/RE closure.

**Numerical integrity:** The reference updates the full firm state and re-evaluates the CES aggregate after every cascade; the Julia 2015 code likewise re-evaluates the CES capital, but does so for *all candidate self-consistent policy vectors*, via ordered breakpoints, then applies ES1. This **existing exact enumerator is more exhaustive** than iterating a single triggered cascade; replacing it would change ES1 and could miss valid equilibria.

## Concrete correspondence

| Component | Reference implementation | Nirei (2015) Julia implementation | Source-based judgement |
|---|---|---|---|
| Parameters and steady aggregate | `setparams.m`, `get_equilibrium.m` bisection, 2-state price thresholds | `calibration.jl` derives closed-form Kbar,Ybar,Xbar and integral A | Different economic models; no transcription error identified |
| Initial heterogeneity | `simulations.m:60–72` conditional inverse-CDF and empirical price normalization | `firms.jl:initialize_firms` draws unif capital gaps `s_i0`, iid initial `a_i0`, fixed `lambda_i` per Appendix C | **Correctly different distribution**; cannot substitute pricing stationary CDF |
| Relative aggregate normalization | sampled price divided by empirical CES P | physical K0 calculated by finite-N CES, sampled individual k_i based on Kbar | **Not directly portable**: nominal P numeraire vs endogenous real physical K; see fixed-point issue below |
| Individual decisions | Calvo reset, upper/lower repricing thresholds | 2015 three-choice `{−1,0,+1}` investment based on `b_i(a_i,t+1)K_t+1` | Threshold *logic* analogous; policy levels not interchangeable |
| Aggregate feedback | reaggregate P after each induced price action | enumerate all self-consistent finite-N CES capital roots via breakpoints | Already handles all reactions; replacing with cascade iteration risks dropping equilibria |
| Equilibrium selection | event direction and trigger determines cascade path | original ES1: minimize `|log(X/Xe)|` across positive-X roots | **Different selection**; cannot import without changing 2015 paper |
| Firm state transition | post-event relative-price re-normalization; wage version adds drift | `solve_equilibrium!` writes next firm capitals and calls `advance_productivity!` only after checks | No proven off-by-one error |
| Aggregate expectations | conditional stationary price-setting thresholds, exogenous inflation in simulations | Appendix A log-linear B\\D and stable eigenvector based on K and lagged mu | The MATLAB code has no 2015 finite-N RE closure to import |
| Additional available code | Hatcher `PF_insert` / `Select_solution`: alternative OBC perfect-foresight regimes and stochastic selection; KhanThomas.jl: HJB + KFE / linearized GE | 2015 ES1 and three-choice policies | Neither is a compatible drop-in algorithm |

## Initial K0 normalization: concrete algebraic obstruction

The original 2015 paper's random initial gap is `s_i0 = log(k_i0/(b_i0 K0))/log(lambda_i)`. The Julia code samples `s_i0` uniformly and forms `k_i0=b_i0 Kbar lambda_i^s_i0`; the exact finite-N CES aggregate returns `K0=Kbar * ( (1/N)sum_i a_i0^(rho/alpha) b_i0^rho lambda_i^(rho*s_i0) )^(1/rho)`.

Thus the *actual* relative gap is `s_i0_actual = s_i0_drawn + log(Kbar/K0)/log(lambda_i)`; it is not literally the sampled uniform gap at finite N. **This is a finite-N consistency discrepancy worth documenting**, but it is **not a proven numerical bug**: imposing exact realized `K0=Kbar` by a common rescaling of all k_i does NOT in general make the originally drawn gaps exactly self-consistent. Simultaneously demanding *independent continuous uniform gaps*, exact finite-N CES, and the identity `k_i=b_i K0 lambda_i^s_i` requires the sampled random sum to equal exactly 1, an event with probability zero. The 2015 numerical paper does not state the necessary finite-N conditional resampling/normalization rule. The 2024 paper's *price numeraire normalization* therefore does not provide a defensible unique correction for 2015 real capital.

## Counterfactual numerical consequences and stop decision

1. **A source-proven 2015 code error was NOT identified.** The audited threshold map and finite-N CES root enumeration are algebraically consistent with the investment policy, and the 2015 timing of iid a_(i,t+1) is consistent with the source. No model-preserving one-line patch is justified by the reference code.
2. **A clear non-portable numerical distinction has been identified:** 2024 code runs *triggered directional avalanches*, whereas 2015's ES1 must compare all finite-N fixed points. This is not an alternative way to guarantee positive aggregate investment. If every fixed point has X<=0, a cascade cannot manufacture a positive fixed point.
3. **A genuinely ambiguous finite-N initial sampling step remains** (physical real aggregate K versus nominal price normalization); do NOT silently impose a new scale constraint, a new gap distribution, or select favorable samples. The paper does not provide a unique model-preserving rule.
4. **Dynamic stochastic RE remains the mathematical bottleneck**, because 2024 MATLAB does not implement the 2015 Appendix A conditional aggregate expectation operator and ES1. The original code's stable-eigenvector continuation is an explicit approximation rather than an independently solved finite-N RE solution.
5. Previous original model: seed 20261101 stops at t5, 20 preregistered seeds: 0/20 at T20, exact CES root oracle PASS. These experiments were NOT rerun. No T1100/Table2 numbers exist.

**No speculative Julia source edit** was made, no new numerical experiment, no Actions workflows, and no new Stage branch. This reference-code comparison does not license claiming the paper has been replicated.

Status of fresh Julia checks: new syntax **not executed** (no Julia code changed); package / new simulation **not executed**. Prior Windows tests on baseline passed. Issue #5 remains OPEN.
