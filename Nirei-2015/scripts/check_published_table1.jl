#!/usr/bin/env julia
# Published Nirei (2015), Table 1, printed p. 972 (PDF page 21).
# This regression MUST fail for the historical erroneous sigma=15.0.
using Nirei2015
p = benchmark_parameters()
@assert p.sigma == 1.5 "Published Table 1 specifies sigma=1.5; current calibration gives $(p.sigma)"
@assert p.alpha == 0.26
@assert p.delta == 0.02
@assert p.beta == 0.99
@assert p.eta == 11.0
@assert p.firms == 110000
c = run_calibration()
@assert c.parameters.sigma == 1.5
println("PASS: published 2015 Table 1 risk aversion sigma=1.5 and calibration")
