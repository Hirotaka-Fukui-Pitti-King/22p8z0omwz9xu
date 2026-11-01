# Nirei (2015), Issue #5 — Stage 11: global equilibrium domain audit

**Status:** Research audit, not a claim that the published 1100-period simulation or Table 2 has been replicated. No production Julia code, calibration, or tests have been changed.

**Primary source:** Makoto Nirei (2015), *An Interaction-Based Foundation of Aggregate Investment Fluctuations*, *Theoretical Economics* 10(3), 953–985, DOI [10.3982/TE1611](https://doi.org/10.3982/TE1611). Page and equation identifiers below refer to the final published article, not earlier working papers.

**Code baseline:** Stage 10 branch \`research/issue5-stage10-fixed-point-feasibility\`, commit \`70c1674b4b141b536d36dbb043ed7dcc30fae197\`; original benchmark N=110000, seed=20261101.

## 1. Determinate findings and limitations

The finite-firm three-action investment model has a **positive aggregate-capital fixed point** for every strictly positive finite firm state and strictly positive future-productivity profile. That statement alone does **not** imply a fixed point with positive aggregate *signed net investment*.

Under the **literal Appendix C** independent Gaussian distribution of log productivity, and conditional on any fixed positive finite predetermined firm state and any finite N, there is a positive-probability set of future-productivity profiles for which **all** aggregate-capital fixed points produce strictly negative aggregate investment. Hence the Section 3.4 log-investment innovation, and an ES1 selector defined through it, lack a real-valued, everywhere-defined stochastic transition under that numerical shock law.

**Crucial source qualification:** This probability theorem cannot be presented as a theorem under *all* the paper's stated assumptions simultaneously. Sections 2.6 and 3.2 separately impose an almost-sure bounded log-productivity-increment condition to justify the one-notch three-action optimal policy, whereas Appendix C specifies nondegenerate independent Gaussian draws in log productivity. The Gaussian specification does **not** satisfy the bounded-increment condition. The paper does not reconcile these assumptions or specify what to do with nonpositive aggregate realized investment. Accordingly, the exact finite-N stochastic **log-shock construction** is not globally justified by the combination of printed specifications. This does *not* prove that all candidate economic equilibria or all 1100-period sample paths are impossible.

The specific observed t=5 failure is *not* identified as an extreme-tail productivity configuration: independent exact enumeration and direct 110000-firm recalculation show three fixed points and negative X for each. The counterexample below proves failure is *possible* under the Appendix C support, not the statistical mechanism causing this particular realization.

## 2. The source assumptions: two substantive inconsistencies

### 2.1 Three-action optimality versus Appendix C shock support

- Section 2.1 (printed pp. 956–957) introduces a bounded-support productivity process and states that the one-notch threshold optimality uses a bounded change in the optimal capital threshold.
- Section 2.6 (printed p. 964) assumes a strict almost-sure bound on \`|log a_{i,t+1} - log a_{i,t}|\` proportional to \`q * alpha * (1/rho - 1) * log lambda_lower\`; the condition is used to rule out optimal capital adjustments of more than one notch. The paper's note on p. 965 explicitly says that broader productivity support can be accommodated by *expanding the choice set* to finitely many notches.
- Section 3.2 (printed p. 968) **explicitly says** idiosyncratic productivity is iid and satisfies the Section 2.6 boundedness condition.
- Appendix C (printed p. 983) instead assigns \`log a_{i,t}\` a normal distribution with nonzero standard deviation 0.05% and draws productivity independently every period. Independent nondegenerate Gaussians have unbounded differences: for every finite \`B>0\`, \`Pr(|log a_{i,t+1}-log a_{i,t}|>B)>0\`.

The two claims cannot be exactly true at once for a nondegenerate iid Gaussian sequence. The numerical specification is the one implemented by \`src/calibration.jl\` and \`src/firms.jl\`. Broadening or truncating the shock process is **not** an authorized source-equivalent fix.

### 2.2 Indivisibility parameter support

Section 3.2 (p. 968) also states \`lambda_i > 1/(1-delta)\`. Appendix C (p. 983) specifies a normal indivisibility parameter with mean 1.028, standard deviation 0.004, truncated at two standard deviations, implying support [1.020, 1.036]. The quarterly delta=0.02 gives \`1/(1-delta)=1.020408163...\`, so a nonempty part of this calibration support violates the strict Section 3.2 restriction.

Earlier controlled Issue #5 diagnostics showed that merely enforcing the strict lower support condition did **not** repair the original t=5 failure. Thus this source inconsistency is real but is *not identified as its unique cause*.

### 2.3 Continuum / finite-N distinction

Section 3.3 (p. 970) says agents forecast with a stationary-uniform-gap **continuum approximation** while finite-N realized aggregate capital and investment are determined from Equations (13) and (14). Appendix B (p. 982) derives Equation (18) by integrating a **one-sided gross investment** response under uniform gaps. The latter is neither a general finite-N identity for realized signed net investment nor a global theorem of exact finite-N conditional rational expectations.

Equation (13) prints a productivity vintage \`a_{i,t}\` even though the adjacent equilibrium construction, Eq. (15), Eq. (16) and Section 3.4 use realized \`a_{i,t+1}\` for the next-period capital decision. This is a genuine printed timing inconsistency, not by itself evidence of a coding error; the implementation follows the unambiguous Section 3.4 protocol.

## 3. Finite-firm equilibrium definitions

Fix \`N<infinity\`, \`0<rho<1\`, \`0<alpha<1\`, \`0<delta<1\`, \`A>0\`, current \`k_i>0\`, finite \`1<lambda_i<infinity\`, and next productivities \`a_i>0\`. Put \`d_i=(1-delta)k_i\` and

\`\`\`text
b_i(a_i) = A * [a_i^(rho/alpha) * (lambda_i^rho-1)/(lambda_i-1)]^(1/(1-rho)).
ell_i = d_i / (lambda_i*b_i)
u_i   = d_i / b_i.
\`\`\`

The one-notch policy is \`kappa_i(K)=-1\` for \`K<=ell_i\`, \`0\` for \`ell_i<K<=u_i\`, and \`+1\` for \`K>u_i\` (the choice exactly at an indifference threshold is inconsequential for continuous productivity draws but must be respected by the enumerator).

Define

\`\`\`text
Phi(K) = [ (1/N) sum_i a_i^(rho/alpha) (d_i lambda_i^kappa_i(K))^rho ]^(1/rho)
X(K)   = (1/N) sum_i d_i (lambda_i^kappa_i(K)-1).
E      = { K>0 : Phi(K)=K }
E+     = { K in E : X(K)>0 }.
\`\`\`

The \`Phi\` and \`X\` response functions are finite, nondecreasing step functions. Their only jumps are upward. \`Phi(0+)>0\` and \`sup_K Phi(K)<infinity\`. Hence \`R(K)=Phi(K)-K\` is positive near zero and negative for large K. Within each breakpoint interval it has slope -1, and at the breakpoints it can jump only upward. Moving from positive to negative therefore requires a **continuous zero crossing** within a breakpoint interval (or at an included endpoint). Thus \`E\` is nonempty and finite for each such positive profile; an appeal to Brouwer continuity would be invalid.

Since \`X(K)\` is nondecreasing, let \`K_X\` denote the first switching boundary after which \`X(K)>0\`. With the actual threshold-tie convention \`X(K_X)<=0\`, and \`X(K)>0\` iff \`K>K_X\`. If \`Kmax=max E\`, then

\`\`\`text
E+ is nonempty   iff   Kmax > K_X.
\`\`\`

This is an exact finite-firm criterion, NOT a Jensen approximation or a test of equilibrium selection.

## 4. Proposition: Appendix C iid lognormal shock implies failure with strictly positive probability

**Proposition.** Fix any deterministic, positive, finite current firm state \`(k_i,lambda_i)_{i=1}^N\`, \`N<infinity\`, and finite parameters in Section 3 with \`0<rho<1\`. Suppose \`a_{i,t+1}\` are conditionally independent, identically distributed, *nondegenerate lognormal* random variables with positive density on all of \`(0,infinity)\`, as in Appendix C. Under the **specified three-action finite-firm reaction system**, conditional on this state,

\`\`\`text
Pr(E+ = empty | state) > 0.
\`\`\`

**Proof.** Write

\`\`\`text
m = rho / [alpha*(1-rho)] > 0,
c_i = A * [(lambda_i^rho-1)/(lambda_i-1)]^(1/(1-rho)) > 0,
b_i(a_i) = c_i * a_i^m.
\`\`\`

Under every possible policy \`kappa_i in {-1,0,1}\`, \`k_{i,t+1}<=d_i*lambda_i\`. For a future-productivity vector whose every component satisfies \`0<a_i<epsilon\`,

\`\`\`text
Phi(K) <= M * epsilon^(1/alpha),
M = [(1/N)sum_i (d_i*lambda_i)^rho]^(1/rho) > 0.
\`\`\`

The lower switching threshold satisfies

\`\`\`text
ell_i = d_i/(lambda_i*c_i*a_i^m)
      >= [d_i/(lambda_i*c_i)] * epsilon^(-m)
      >= L * epsilon^(-m),
L = min_i d_i/(lambda_i*c_i) > 0.
\`\`\`

Choose \`epsilon>0\` so small that

\`\`\`text
M * epsilon^(1/alpha) < L * epsilon^(-m).
\`\`\`

Every capital fixed point then obeys \`K=Phi(K)<min_i ell_i\`. Every firm consequently chooses \`kappa_i(K)=-1\`. Its aggregate signed net investment equals

\`\`\`text
X(K) = (1/N)sum_i d_i*(lambda_i^(-1)-1) < 0.
\`\`\`

The event \`{0<a_i<epsilon for every i}\` has positive conditional probability \`[F_a(epsilon)]^N>0\` for finite N under independent nondegenerate lognormals. It is contained in \`{E+=empty}\`. This proves the claim. QED.

**Scope.** The proof is exact for the finite-N *specified one-notch reaction system plus Appendix C shock law*. It does not establish the probability under the *separate and incompatible* Section 2.6 a.s.-bounded-increment assumption, does not show the extremely small all-low-a event is quantitatively important, and does not show that the 3-choice rule remains a globally optimal dynamic firm policy outside that condition. It is a domain **counterexample**, not a claimed proof of nonexistence of a full general equilibrium for all possible extensions.

## 5. The log-investment innovation has only a partial domain

Eq. (14), p. 970, defines signed net firm investment relative to depreciated capital: \`X_t=(1/N)sum_i (k_{i,t+1}-(1-delta)k_{i,t})\`. It permits divestment, so \`X_t<0\` is arithmetically possible. Eq. (12) is \`Y_t=C_t+X_t\`; negative X does not **by itself** contradict the resource identity, which must additionally satisfy positive consumption and other equilibrium conditions.

Eq. (19), p. 970, writes \`X_t = E_{t-1}[X_t]*exp(epsilon_t)\`, and Section 3.4, p. 971, defines epsilon as a log difference of realized and expected investment. A real-valued log shock requires \`X_t>0\` and \`E_{t-1}[X_t]>0\`. The implemented Appendix A forecast is positive by level exponentiation. The paper's Section 3.4 ES1 minimizes the absolute log investment surprise among feasible fixed points, but provides **no fallback** when every fixed point has \`X<=0\`.

Thus the simulation's specified transition from state, future productivity and an ES1 choice to a real log shock is a **partial map**, defined on

\`\`\`text
D = { (state,a_next) : E+(state,a_next) is nonempty
                         and expected X(state)>0 }.
\`\`\`

The Appendix C conditional transition law assigns \`D^c\` positive probability at every fixed positive finite state. In particular the model as *stated and implemented* does not supply a real-valued, a.s.-defined stochastic process for \`epsilon_t\`. An extension that selects some capital fixed point with \`X<=0\` can be defined mathematically but is **not** the paper's log-ES1 shock system, and the original expectation and endogenous continuation would then need to be re-solved.

A random trajectory may remain in \`D\` for 1100 periods for a suitable seed; that observation cannot establish a globally a.s.-defined transition. Nor does one failed seed falsify every reported simulation. Exogenous shock rejection instead conditions on survival/feasibility and changes the transition kernel.

## 6. Appendix A and exact rational expectations

Appendix A Proposition 5 (printed pp. 980–982) log-linearizes the continuum expectation system, stacks \`B E_{t-1} z_t=D z_{t-1}\`, and proves a local saddle path under the investment/output restriction \`a_x<alpha\` (proposition statement presents a weak inequality). It also writes a first-order innovation update \`tilde K_0=E_{-1}tilde K_0+delta*epsilon_0\`. The following are distinct:

1. A deterministic steady state / continuum uniform-gap *forecast approximation*.
2. A locally determinate linear rational-expectations system conditional on its valid first-order log innovation.
3. The actual finite-N equilibrium capital correspondence under Eq. (13).
4. A globally well-defined *nonlinear finite-N stochastic rational-expectations equilibrium* with a total transition rule and conditional expectations over every permitted future shock.

Proposition 5 shows (2), not (4). In \`src/expectations.jl\` the code constructs \`F=B\\D\` and its stable eigenvector \`h\`. In \`src/simulation.jl\`, realized marginal utility is **imposed** via \`tilde mu_t=h*tilde K_{t+1}\`, and contemporaneous macro levels are exponentiated first-order log solutions. A stable-manifold relation may be justified in an appropriate *stochastic linear RE solution*, if innovations lie on its stable subspace and relevant conditional expectations satisfy the model restrictions; it does not by itself prove finite-N nonlinear Euler/resource/expectation consistency.

An unconditional \`E[X_t^{ES1}|state]\` cannot be inferred from MC conditioned on \`E+!=empty\` unless a well-defined choice on \`E+=empty\` is supplied. An alternative signed-X capital equilibrium selection could define an unconditional signed-X random variable, but that would be a **different choice rule** from Section 3.4 log-ES1. Mean-zero log investment shock should not be imposed mechanically: \`E[log X]!=log E[X]\` in general.

## 7. Existing exact numerical counterexample and its implications

**Verified on user's Windows Julia 1.12.6** and in [Stage 10 GitHub Actions 41259640042](https://github.com/Hirotaka-Fukui-Pitti-King/22p8z0omwz9xu/actions/runs/41259640042) at commit \`70c1674\`: syntax 18 Julia files, package import, independent N=48 fixed-point oracle across 4 seeds, benchmark direct reaggregation and 32-draw feasibility diagnostic **PASS**. Overall workflow **FAILS** on original unchanged no-rejection 20-period benchmark at t=5; do not confuse this with a diagnostic failure.

**Original t5 profile** \`N=110000\`, seed \`20261101\`:

| ascending fixed capital | realized aggregate net X |
|---:|---:|
| 7.922286680029842 | -0.05758656336226558 |
| 7.922473044264564 | -0.0574001973735959 |
| 7.922475106160655 | -0.057398135458978324 |

Original t5 \`K_expected_Eq18=8.146420155827748\`, \`X_fixed(K_expected)=0.15990391297028964\`. All roots were independently verified (max capital residual \`2.842170943040401e-14\`, investment residual \`4.85722573273506e-16\`).

The independently swept strict positive-X threshold \`K_X=7.980278122719776\`, whereas \`Kmax=7.922475106160655\`, so \`Kmax-K_X=-0.05780301655912101\`. Thus no admissible positive-X fixed point exists. The 32-draw **separate conditional iid** diagnostic from the *original frozen t5 firm state* yields 31 feasible, one infeasible, \`p_hat=0.96875\`, 95% Wilson interval \`[0.8425573617998744,0.9944621398359969]\`. This is **not** an unconditional full-path failure probability or a precise rare-event estimate.

The initial investment forecast wedge decomposition previously showed ~94% of the conditional shortfall between fixed forecast-capital firm response and imposing the finite-N equilibrium. ES1 chose the largest positive-investment root in every feasible tested initial state, excluding the claim that simply choosing a higher positive-X root would repair the benchmark failure.

## 8. Source-to-code reconciliation and decision

| Object | Published source | Current implementation | Assessment |
|---|---|---|---|
| Signed aggregate net investment | (14), p. 970 | \`src/equilibrium.jl\` signed \`d_i*(lambda_i^kappa-1)\` | Consistent |
| Three-choice threshold | (10), (16), pp. 969–970 | \`src/firms.jl\`, \`src/equilibrium.jl\` | Formula consistent; global optimality requires the bounded-increment assumption absent from Appendix C shocks |
| Productivity law | Section 2.6/3.2 vs Appendix C | \`src/firms.jl\` iid lognormal | Matches Appendix C; **not** Section 2.6 almost-sure bounded increments |
| Lambda distribution | Section 3.2 vs Appendix C | \`src/calibration.jl\` truncated Gaussian [1.020,1.036] | Matches Appendix C; violates strict Section 3.2 support on [1.020,1.020408...] |
| CES capital | (13), p. 970 | \`src/firms.jl\`, \`src/equilibrium.jl\` | Correct power aggregation, but printed \`a_{i,t}\` vs context \`a_{i,t+1}\` is an unresolved notation inconsistency; implementation follows Section 3.4 timing |
| Log investment innovation/ES1 | (19), Section 3.4, pp. 970–971 | \`src/equilibrium.jl\` positive-X log objective and explicit domain failure | Correct domain enforcement; **no source-defined fallback** |
| Approximate expected capital | (18), Appendix B | \`src/expectations.jl\`, \`src/simulation.jl\` | Local continuum/first-order expectation, not an exact finite-N expectation identity |
| Realized marginal-utility continuation | Appendix A Proposition 5 | \`src/simulation.jl\` uses \`tilde mu=h*tilde Knext\` | Imposed linear stable continuation, not independently established global nonlinear RE |

**Decision:** Outcome **C** (only local/conditional log-ES1 simulation specified by these equations), with an unresolved original-implementation-information caveat under **D**. No new *source-justified* production correction or missing explicit global positive-X theorem has been located in the published paper or Appendix C. The exact original author code / special fallback is unknown; **do not invent it or claim it did not exist**.

Possible *extensions*, NOT replication repairs, include a shock-support/state-invariance restriction guaranteeing \`E+!=empty\`, a total capital equilibrium selection with signed-X dynamics and a newly solved household continuation, or a redefined investment index/innovation. Each changes the economic object; none is authorized as an invisible bug fix.

**No Julia code changed in Stage 11.** Previous Windows tests remain valid for unchanged code; commit of this mathematical audit need not be described as new local numerical execution. Keep Issue #5 **OPEN**; do not run Issue #6 or claim 1100-period non-rejection replication.
