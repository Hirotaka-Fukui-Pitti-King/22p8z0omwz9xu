# Stage 2 author-source audit, timing and verification

## Sources reviewed

The original public author code (all five full source files were consulted):
- https://github.com/kurtmitman/HANKSOME/blob/main/Main.m
- https://github.com/kurtmitman/HANKSOME/blob/main/MIT_transition.m
- https://github.com/kurtmitman/HANKSOME/blob/main/Transition/solve_transition_bis.m
- https://github.com/kurtmitman/HANKSOME/blob/main/Transition/backsolve_egm.m
- https://github.com/kurtmitman/HANKSOME/blob/main/Figures.m

Julia resources: https://julialang.org/ and https://docs.julialang.org/en/v1/manual/performance-tips/ . The LinearAlgebra standard-library UUID is 37e2e46d-f89d-539d-b4ee-838fcccc9c8e; verified from Julia-language packages using the standard library.

## Source to Julia correspondence

- Main.m: 200 quarterly credit shocks B[t] recursively decrease by 0.074 Y_initial 0.98^(t−1). Initial B=0; terminal B corresponds to shock[200].
- solve_transition_bis.m: K[1]=initial.K; unknown K[2:T] and L[1:T] form 399 unknowns. Initial q0=1, ending q[T]=1; investment and q calculated from K path. The authors' period-t r_t uses current rK and current q, divided by previous q.
- solve_transition_bis.m: period income includes home-firm profits Y/eps_f and investment-firm profits (q−1)*I−Psi. Initial aggregate wealth equals q0 K0+B0; interior period-end aggregate wealth equals q[t] K[t+1]+B[t]; terminal wealth equals K_final+B[T].
- solve_transition_bis.m: period terms-of-trade condition equalizes budget-implied home consumption and Y−I−Psi−D*x^(−theta_star). Given theta=1, this scalar equation is solved algebraically for x, but no economic restriction is omitted.
- backsolve_egm.m: terminal coarse and fine consumption policies are final steady state policies. Backward current consumption uses beta*(1+r[t+1]) times next period expected marginal utility and current budget gross return 1+r[t]. PCHIP is applied to endogenous consumption and current-asset nodes directly. Borrowing-constrained saving is at the borrowing limit.
- backsolve_egm.m: initial joint household distribution is the initial SS distribution; forward simulation splits next-asset choice between adjacent asset knots and applies the productivity Markov transition matrix. Consumption and savings are aggregated using the contemporaneous household distribution.
- solve_transition_bis.m: residuals are capital-market clearing for t=1,...,T-1 and union labor clearing for t=1,...,T.
- Figures.m: Figure 1 includes BOTH initial credit expansion and later surprise contraction; Stage 2 is only the first component. Figure 2 displays home and foreign consumption, output and investment. Figure 3 displays trade prices, real wage, nominal exchange rate and real interest rate, compared with initial steady state.

## Important unresolved author-code discrepancy

Main.m sets varphi=0.5 and the steady-state union FOC uses L^varphi. But solve_transition_bis.m calculates the implied labor supply without the 1/varphi power. At an initial stationary equilibrium with L not equal to one, author-code stationary-limit labor residual can therefore be nonzero. This discrepancy must be reported, not masked by modifying the published transition solver while claiming author-code replication. The smoke script reports the initial-vs-implied labor levels.

## Solver and acceptance

Julia's full root solver uses a dense finite-difference Jacobian, damped updates, rank-one Broyden corrections and periodic recomputation. It is not the same numerical optimizer as MATLAB's fsolve, but it solves the same nonlinear GE equation system. All residual evaluations recompute the backward household policies and forward household distribution. Invalid numerical trial states are rejected; no state is clipped merely to force convergence.

Independent diagnostics measure asset, labor, goods, capital accumulation, household budget, distribution mass, dynamic Euler/KKT and terminal K/L errors. The goods diagnostic uses *actual household consumption* from the forward distribution, not the budget quantity used to solve the trade price. The source Figure 1 full two-path comparison cannot be reproduced in Stage 2, because the later surprise contraction is outside scope; the SVG title explicitly states expansion only.

## Status

The Stage 2 Julia source, tests and execution/figure scripts are committed but have **not** been run in Julia on the user's Windows PC. Syntax not verified; Julia import not verified; numerical convergence not verified; figures not yet generated; comparison with author numerical paths not verified. Do not claim complete Figure 1 replication.


### 2026-11-02: Windows Stage 2 first-run failure and targeted correction

Owner Windows 11 / Julia 1.12.6 result: scripts/check_syntax.jl printed PARSE OK on 26 files, but its previous parser did not reject returned parse-error expressions. test/runtests.jl passed 2,420 Stage 1 tests, then Stage 2 test parsing failed at test/test_stage2.jl:53 because view(...,:,:,end) uses end outside indexing syntax.

The eight-period smoke completed backward EGM and forward distribution, including mass error 1.4432899320127035e-14, but its stationary candidate was NOT an asset-clearing GE (beta was imposed 0.9832). Its 0.009541840016043829 asset and 0.0012405672215367591 labor residuals are therefore not proof of a bad dynamic model.

The 200-period simulation finished three stationary GE solves but FAILED before the first root iteration: EGM saving at a far-right fine asset node was 29759.54040678395, above the upper distribution grid. Unlike the original author's backsolve_egm.m, the Julia backward household block erroneously imposed the distribution upper asset bound on all fine-grid policy nodes, including nodes with zero population mass.

The correction follows author EGM: every consumption stays strictly positive and savings satisfies the lower borrowing bound. No upper-grid bound is placed on policy evaluations at unoccupied high-wealth nodes. The forward distribution applies the upper-grid condition to all **positive-mass** states and throws an error for any true out-of-grid decision. Nothing is silently clipped or projected outside the distribution grid. A direct regression covers both occupied and unoccupied off-grid choices.

The parser now rejects error/incomplete AST nodes, the test replaces invalid end by explicit size(...,3), and the stationary-limit smoke calibrates a true stationary GE before constructing a constant-B path. The known author labor exponent inconsistency is preserved, not silently altered.

PATCH STATUS: committed but Julia syntax, full tests, eight-quarter true-stationary test, 399-variable nonlinear convergence and generated SVGs remain UNVERIFIED on the owner's PC.


### 2026-11-02: Windows 200-quarter nonlinear solve succeeded; terminal goods check isolated

Windows 11 / Julia 1.12.6 **owner-provided local execution log** for the previous commit:

- 26 Julia files passed improved Meta.parseall syntax checks.
- 2,420/2,420 stationary tests and 30/30 Stage 2 transition tests passed.
- True stationary GE constant-B eight-quarter limit: max scaled GE residual 1.2327219356400576e-7, max asset 6.362651475204854e-8, max labor 1.2327219356400576e-7, distribution mass deviation 7.327471962526033e-15.
- Full T=200, 399-variable nonlinear root **converged** in 7 iterations with 409 residual evaluations, max scaled 399-equation residual 3.642975265893966e-7.
- Independent asset market maximum level gap 9.364358355412605e-7; labor maximum 3.642975265893966e-7; capital accumulation gap 0; maximum household budget gap 1.3304388854170658e-16; maximum distribution-mass error 1.2212453270876722e-14.
- Existing independent goods-market maximum reported **0.2633831568791398**, but it pooled **all t=1:200** rather than only the interior dates for which individual household savings and consumption are part of the original author's transition closure.
- Maximum dynamic Euler/KKT violation 0.005496693743120638 (mass weighted mean 1.9816793256669015e-5). This does not meet the former strict 1e-6 maximum accuracy target.
- Terminal capital discrepancy abs(K[200]-K_final)=0.42899657754234255, terminal labor discrepancy 0.0071350035667030465. Both must be explicitly reported; they are finite-horizon approximations, not hidden by the root solver.
- CSV and SVG files were emitted before the original script discovered a failing goods diagnostic. Those files are **unaccepted provisional outputs**, not verified replications.

**Confirmed source-level endpoint mismatch:** the original author `backsolve_egm.m` aggregates household consumption and saving at dates 1:TT-1 and the original `solve_transition_bis.m` sets `c_t=[c_t';c_t(end)]`, reusing consumption of TT-1 for the last labor FOC. The T-th consumption policy and distribution are terminal evaluation objects, not an additional independent market-clearing equation. The Julia diagnostic mistakenly included T-th household consumption in the maximum independent goods-market residual, imposing a condition beyond the authors' 399 equations.

**Correction committed:** preserve the independent household-consumption goods check for **interior dates 1:199 only**; separately output `terminal_goods_gap_unconstrained`, `terminal_home_budget_gap_unconstrained`, and `terminal_wealth_gap_unconstrained`. Record `max_goods_market_period` and `max_interior_home_budget_gap` to determine if the 0.263383 gap was all at the endpoint. Do NOT loosen any tolerance, adjust the labor FOC, or change household EGM. Generate new CSV/SVG only after the independent interior GE acceptance gate passes. Regression test deliberately corrupts the T-th reported aggregate consumption and confirms the interior goods residual is unchanged, while the terminal diagnostic responds.

**Verification status:** this newest diagnostic patch is not yet syntax-checked or executed on the owner's Windows PC. It is still unknown whether the maximum interior goods residual meets 1e-4; the new code will report its period. The T=200 root convergence reported above is established for the previous commit, not the patched diagnostic.


### 2026-11-02: full Stage 2 independent aggregate GE validation PASSED on owner PC

**Evidence**: Owner-provided Windows 11 / Julia 1.12.6 PowerShell output from commit `15a443508eddadcc30512b89310d74e006d3e061`.

| Item | Observed PC result |
|---|---:|
| Recursive Julia syntax check | 26/26 passed |
| Stage 1 tests | 2420/2420 passed |
| Stage 2 tests | 37/37 passed |
| 8-quarter constant-B stationary limit, scaled GE residual | 1.2327219356400576e-7 |
| 8-quarter interior goods-market discrepancy | 1.797111210954938e-7 |
| Nonlinear unknowns / transition periods | 399 / 200 |
| Nonlinear iterations / function evaluations | 7 / 409 |
| Maximum normalized 399-equation residual | 3.642975265893966e-7 |
| Maximum independent asset-market level discrepancy | 9.364358355412605e-7 |
| Maximum independent labor-market level discrepancy | 3.642975265893966e-7 |
| Maximum interior (t=1,...,199) goods-market discrepancy | 5.52355026384177e-7 (t=1) |
| Maximum capital accumulation discrepancy | 0 |
| Maximum household budget discrepancy | 1.3304388854170658e-16 |
| Maximum probability mass discrepancy | 1.2212453270876722e-14 |
| Maximum dynamic Euler/KKT violation | 0.005496693743120638 |
| Probability-weighted mean dynamic Euler/KKT violation | 1.9816793256669015e-5 |
| Terminal capital discrepancy | 0.42899657754234255 (approximately 1.33% of final K) |
| Terminal labor discrepancy | 0.0071350035667030465 (approximately 0.70% of final L) |
| Terminal goods discrepancy (not an additional market equation) | -0.2633831568791398 |
| Terminal home-budget discrepancy (not an additional market equation) | -0.2633831568791398 |
| Terminal wealth discrepancy (not an additional market equation) | -0.4294718146260834 |

**Decision**: The full 200-quarter, 399-variable nonlinear transition root and independent **interior** asset, labor, and goods-market checks have passed the implemented acceptance criteria. The earlier maximum goods discrepancy of 0.263383 was entirely in the T=200 terminal diagnostic, not the interior markets: max interior goods discrepancy is only 5.52e-7. Terminal conditions remain approximation diagnostics; they must not be described as exactly reached. The dynamic maximum Euler residual exceeds the earlier 1e-6 policy target and remains an accuracy limitation.

Owner's PC successfully generated `results/stage2_transition_expansion.csv`, `results/stage2_equilibrium_diagnostics.csv`, `results/stage2_solver_log.txt` and the three SVG files. **Their new bytes are on the owner's local PC and are NOT yet in the GitHub branch.** Previous GitHub `results/stage2_solver_log.txt` is a stale failure log until explicitly replaced.

### Original Figure 1–3 source mapping audit

Authors' `Figures.m` (public GitHub source, lines 21–195) was re-read after the accepted Windows run:

- Figure 1 plots `-0.25*shock_b_agg/aggregates_initial.y` for the expansion plus an unanticipated contraction path. Current Stage 2 SVG correctly limits its scope to the expansion component; it is **not** a complete Figure 1 replication.
- Figure 2 plots home and foreign consumption, output as `100*(log(variable)-log(initial))`; investment uses `100*(I-Iinitial)/Iinitial`. The current SVG plot formulas follow these transformations.
- Figure 3 plots terms of trade, real wage, inverse terms of trade as the nominal exchange rate in log percentage deviations, and the quarterly real interest rate in percentage **levels** with an initial-SS point. The current SVG plot formulas follow these definitions.
- The author source fixes panel axis ranges and adds vertical markers at the shock date. The Julia SVG currently uses automatic axis ranges, so visual presentation is not an exact typographical copy.
- The author's root repository includes `transition_start.mat` but no committed `solution_figures.mat`. This inspection does not establish quantitative equality of Julia trajectories to the author's MATLAB trajectories. Generated SVG content has not yet been visually inspected by the assistant (local files not uploaded).

**Stage 2 aggregate-equilibrium implementation: locally validated. Figure 1 complete replication: not done (later unexpected contraction is outside Stage 2). Figure 2–3 source variable/normalization correspondence: audited. Quantitative author-series agreement: not verified. Full household Euler precision target: not met.**


### Actual Stage 2 artifacts inspected on GitHub (after owner push)

The owner pushed the successful Windows-run result files to the branch, which were incorporated through merge commit `948d88cd85fdb9f2113f40e0eb8852701b67ddf9`. This **supersedes the earlier claim** that the updated CSV/SVG files were only on the user's PC. Verified public branch files now include:

- `results/stage2_transition_expansion.csv` (200 quarterly observations; 23 columns)
- `results/stage2_equilibrium_diagnostics.csv`
- `results/stage2_solver_log.txt` (updated successful solver log)
- `results/figure_1.svg`, `results/figure_2.svg`, `results/figure_3.svg`.

Independent CSV inspection confirmed the logged market residuals, as well as the source-derived Figure 2–3 transformations. Initial stationary values from `results/stage1_steady_states.csv` were used as denominators. The **first simulated quarter** (CSV quarter=0; the plotted image places it at horizontal tick=1 after the initial steady-state anchor) has:

| Series, original figure's definition | First simulated quarter |
|---|---:|
| Home-good consumption, 100 × log ratio to initial SS | +4.96294% |
| Foreign-good consumption, 100 × log ratio to initial SS | +11.72675% |
| Output, 100 × log ratio to initial SS | -2.09740% |
| Investment, 100 × level percentage difference | +7.29389% |
| Terms of trade, 100 × log ratio | +6.76381% |
| Real wage, 100 × log ratio | +1.03305% |
| Nominal exchange rate, 100 × log ratio | -6.76381% |
| Real interest rate, quarterly percentage level | 3.49331% |

At CSV quarter=99 (outside the first 100 *plotted* points because the image inserts the initial stationary point), domestic consumption is +0.48385%, foreign consumption -1.76878%, output +4.74383%, investment +12.39952% relative to initial SS; terms of trade -2.25263%; wage +2.21410%; nominal exchange +2.25263%; interest 0.95469% per quarter.

**Terminal caution from the actual CSV**: between quarters 198 and 199, terms of trade jump from 0.9821022495 to 0.9131088404, whereas capital remains near 31.837. This is a T-th boundary approximation associated with the reported terminal goods and wealth gaps, not a regular interior market-clearing outcome. The first 100 plotting points do not include this terminal observation.

Original `Figures.m` panel definitions and scaling agree with the plotted Julia formulas. The original Figure 1 includes a later surprise contraction, which Stage 2 does not implement. The original `Figures.m` uses fixed y-axis ranges while the Julia SVG uses dynamic ranges. The **exact original numerical response paths are not supplied in the author's public repository** (the main branch has `transition_start.mat` but not `solution_figures.mat`), so matching the source panel variables is **not** confirmation of numerical equality to the paper. The SVG files were fetched and verified as files/XML, but no graphical rendering inspection was completed in the assistant environment. Keep quantitative paper-figure replication marked as unverified.
