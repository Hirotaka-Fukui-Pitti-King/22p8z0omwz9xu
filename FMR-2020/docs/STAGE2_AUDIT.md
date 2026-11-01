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
