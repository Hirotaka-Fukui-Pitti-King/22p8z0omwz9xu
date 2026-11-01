#!/usr/bin/env julia
# Equations (7),(8),(11),(12) exact-level regression. NOT conditional RE.
using Nirei2015

# Function scope avoids Julia's global soft-scope ambiguity in try/catch.
function main()
cal=run_calibration()
s=cal.steady_state
p=cal.parameters
@assert p.sigma==1.5
base=exact_static_levels(cal,s.K,s.X,s.mu)
for key in (:Y,:C,:L,:w,:mu,:R)
    actual=getproperty(base,key)
    expected=getproperty(s,key)
    @assert isapprox(actual,expected;rtol=5e-10,atol=5e-12) (key,actual,expected)
end
off=exact_static_levels(cal,s.K,1.08*s.X,s.mu)
for key in (:resources_residual,:labor_residual,:production_residual,:pricing_kernel_residual)
    @assert abs(getproperty(off,key))<1e-10 key
end
@assert off.C>0.0 && 1.0-p.psi*off.L^p.zeta>0.0
limit=s.K^p.alpha*p.psi^(-(1.0-p.alpha)/p.zeta)
was_domain_error=false
try
    exact_static_levels(cal,s.K,1.01*limit,s.mu)
catch err
    was_domain_error=err isa DomainError
end
@assert was_domain_error
println("PASS: exact nonlinear steady-state/off-steady level conditions and infeasible-X domain")
println("NOT TESTED: nonlinear conditional pricing expectations or full stochastic RE")
end

main()
