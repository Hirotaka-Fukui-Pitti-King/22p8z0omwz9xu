
using Test
using FMR2020

@testset "Stage 3 shock and balance-sheet mechanisms" begin
    @testset "Unexpected reversal source timing" begin
        B=collect(1.0:1.0:200.0)
        result=credit_reversal_path(B,-4.0;shock_start=41,periods=200)
        @test length(result)==200
        @test result[1]==B[40]
        @test result[end]==-4.0
        @test isapprox(result[2],-4+(40+4)*(1-1/199)^4;atol=1e-12)
        @test_throws ArgumentError credit_reversal_path(B,-4.0;shock_start=1)
    end
    @testset "Positive trade-price roots with same-sign far endpoints" begin
        # This is the actual coefficient vector implied by the failed
        # T=200 initial aggregate guess in the owner's Windows log.
        C=-2.7782230139198565
        E=3.8257097540368417
        D=1.0149463596999195
        theta_star=3.0
        f(x)=C+E/x-D/x^theta_star
        xpeak=sqrt(3*D/E)
        @test f(0.001)<0
        @test f(1000.0)<0
        @test f(xpeak)>0
        root=FMR2020.foreign_debt_trade_root(C,E,D,theta_star)
        @test isapprox(root,0.7865895839646455;atol=1e-9)
        @test root<xpeak
        @test abs(f(root))<1e-9
        roots=FMR2020.foreign_debt_trade_roots(C,E,D,theta_star)
        @test isapprox(roots.low,root;atol=1e-11)
        @test roots.high !== nothing
        @test isapprox(roots.high,1.0379166453624395;atol=1e-9)
        @test roots.low<roots.peak<roots.high
        @test abs(f(roots.high))<1e-9
        # Opposite signs at two arbitrary far endpoints are NOT
        # required, but genuine lack of any positive root still fails.
        @test_throws DomainError FMR2020.foreign_debt_trade_root(
            -10.0,1.0,1.0,3.0)
    end

    @testset "Regularized nonlinear direction remains descent near singular Jacobian" begin
        # The user's full GE stopped because the 399-variable Broyden
        # Jacobian was not refreshed after many near-zero damped steps.
        # The LM rescue direction must lower the local least-squares
        # objective even for a nearly singular finite-difference J.
        J=[1.0 1.0; 1.0 1.0+1e-11]
        residual=[0.02,-0.01]
        for λ in (1e-6,1e-4,1e-2,1.0,100.0)
            step=FMR2020.stage3_lm_direction(J,residual,λ)
            @test all(isfinite,step)
            grad=transpose(J)*residual
            @test sum(grad.*step)<0
        end
        @test_throws ArgumentError FMR2020.stage3_lm_direction(
            J,residual,0.0)
        @test_throws DimensionMismatch FMR2020.stage3_lm_direction(
            ones(2,3),residual,1e-3)
    end

    @testset "Zero exchange-rate shock preserves household net wealth" begin
        grid=Float64[-2,-1,0,1,2,3,4]
        mu=zeros(7,2)
        mu[3,1]=0.1
        mu[4,1]=0.3
        mu[5,2]=0.6
        dest=similar(mu)
        info=revalue_foreign_debt!(dest,mu,grid,2.0,1.0,-0.5,
                                    1.2,1.2,1/16)
        @test abs(info.mass_after-1)<1e-12
        @test abs(info.net_wealth_before-info.net_wealth_after)<1e-12
        @test maximum(abs.(dest.-mu))<1e-12
        @test isapprox(info.capital+info.bonds,1.5;atol=1e-12)
        @test_throws ArgumentError revalue_foreign_debt!(
            mu,mu,grid,2.0,1.0,-0.5,1.2,1.2,1/16)
    end
    @testset "Exchange-rate change revalues liabilities and conserves mass" begin
        grid=Float64[-3,-2,-1,0,1,2,3,4]
        mu=zeros(8,1)
        mu[3,1]=0.5  # borrower with negative foreign bond holdings
        mu[6,1]=0.5
        dest=similar(mu)
        info=revalue_foreign_debt!(dest,mu,grid,2.0,1.0,0.0,
                                    1.0,2.0,0.0)
        @test isapprox(sum(dest),1.0;atol=1e-12)
        @test isapprox(info.net_wealth_after,
                       info.capital+info.bonds/2;atol=1e-12)
        @test minimum(dest)>=0
        @test isapprox(info.net_wealth_after-info.net_wealth_before,
                       0.25;atol=1e-12)
    end
    @testset "Positive probability leaving full asset grid errors" begin
        grid=Float64[-1,0,1]
        mu=zeros(3,1)
        mu[1,1]=1.0
        dest=similar(mu)
        @test_throws DomainError revalue_foreign_debt!(
            dest,mu,grid,1.0,1.0,0.0,1.0,0.1,0.0)
    end
end
