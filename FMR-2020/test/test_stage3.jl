
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
