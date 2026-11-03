using Test
using FMR2020

include(joinpath(@__DIR__,"..","scripts","run_economic_replication.jl"))

@testset "Saved economic market-root resumption is checked, not trusted" begin
    p=default_parameters(n_coarse=61,n_fine=241,
                         egm_tol=1e-9,dist_tol=1e-9)
    steady=solve_three_steady_states(p;verbose=false)
    ctx=stage3_expansion_context(steady)
    n=length(ctx.economy.B)
    @test n==200
    K=fill(ctx.economy.initial.K,n)
    L=fill(ctx.economy.initial.L,n)
    terms=fill(ctx.economy.initial.terms,n)
    mktempdir() do d
        csv=joinpath(d,"synthetic.csv")
        function export_path(;wrong_B=false)
            open(csv,"w") do io
                println(io,
                    "period,B,K,L,Y,I,q,r,w,terms,foreign_consumption,home_consumption,wealth,nominal_exchange_rate")
                for t in 1:n
                    B=ctx.economy.B[t]+((wrong_B && t==2) ? 0.25 : 0.0)
                    println(io,join((t,B,K[t],L[t],1.0,0.1,1.0,0.01,0.5,
                                     terms[t],1.0,1.0,1.0,1.0),","))
                end
            end
        end
        export_path()
        v=read_economic_root_guess(csv,ctx)
        @test length(v)==3n-1
        @test v[1:n-1]==K[2:end]
        @test v[n:2n-1]==L
        @test v[2n:3n-1]==terms
        fake=(K=copy(K),L=copy(L),Y=ones(n),I=fill(0.1,n),
              wage=fill(0.5,n),rate=fill(0.01,n),
              terms=copy(terms),aggregate_home=ones(n))
        same=compare_boom_horizons(csv,fake)
        @test same.maximum==0.0
        @test same.first_periods==40
        fake.K[5]+=0.03
        different=compare_boom_horizons(csv,fake)
        @test isapprox(different.by_variable["K"],
                       0.03/max(1.0,abs(K[5]));atol=1e-12)
        export_path(wrong_B=true)
        @test_throws ArgumentError read_economic_root_guess(csv,ctx)
        rm(csv)
        @test read_economic_root_guess(csv,ctx)===nothing
    end
end
