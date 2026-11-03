using Test
using FMR2020

include(joinpath(@__DIR__,"..","scripts",
                 "stage8_depreciation_continuation.jl"))

@testset "Stage 8 source numerical guess/metadata integrity" begin
    input=normpath(joinpath(@__DIR__,"..","results","author_reference"))
    capital_labor=source_vector(input,"leveraged_candidate")
    terms=source_vector(input,"leveraged_source_terms")
    inherited=source_vector(input,"shock_inherited_metadata")
    @test length(capital_labor)==399
    @test length(terms)==200
    @test length(inherited)==7
    @test all(isfinite,capital_labor)
    @test all(isfinite,terms)
    @test all(>(0),terms)
    @test terms[1]<inherited[4]
    @test isapprox(terms[1],0.90985851141107077;atol=1e-12)
    @test isapprox(inherited[4],1.0258494514359664;atol=1e-12)
    @test length(vcat(capital_labor,terms))==599
end

@testset "Stage 8 rejected-candidate state is serializable without false convergence" begin
    mktempdir() do d
        name=joinpath(d,"trial.csv")
        rows=Pair{String,Any}[
            "source_guess_only"=>true,
            "equilibrium_root_status"=>"NOT CERTIFIED",
            "failure_type"=>"PositiveMassEconomicInfeasibility",
            "offgrid_wealth"=>-0.25,
        ]
        stage8_store_metrics(name,rows)
        @test isfile(name)
        contents=read(name,String)
        @test occursin("equilibrium_root_status,NOT CERTIFIED",contents)
        @test occursin("failure_type,PositiveMassEconomicInfeasibility",contents)
        @test !occursin("CONVERGED",contents)
    end
end

println("STAGE8 LIGHT UNIT TESTS PASS (no GE residual/Jacobian/solve)")
