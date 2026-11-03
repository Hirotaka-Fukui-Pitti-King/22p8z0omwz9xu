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

@testset "Stage 8 restart reads saved 599-vector and checks residual provenance" begin
    n=599
    pre=1.0271423614420732
    x=ones(Float64,n)
    x[1]=29.1
    x[200]=1.02
    x[400]=0.93
    residual=fill(0.0015,n)
    mktempdir() do d
        path=joinpath(d,"source_initialized_UNCONVERGED.checkpoint")
        FMR2020.stage3_write_checkpoint(path,x,residual,5)
        state=stage8_read_checkpoint(path,n)
        @test state.iterations==5
        @test state.x==x
        @test state.stored_max==0.0015
        @test stage8_verify_resume_residual(state.stored_max,0.0015)
        @test_throws ErrorException stage8_verify_resume_residual(
            state.stored_max,0.0025)
        @test_throws ArgumentError stage8_read_checkpoint(
            joinpath(d,"absent.checkpoint"),n)
        @test_throws DimensionMismatch stage8_read_checkpoint(path,n+1)
        trace=joinpath(d,"accepted_full599_trace.csv")
        stage8_append_accepted_step(trace,"--resume",5,1,x,residual,pre,200)
        @test isfile(trace)
        lines=readlines(trace)
        @test length(lines)==2
        headers=split(lines[1],",")
        fields=split(lines[2],",")
        @test length(headers)==length(fields)
        @test fields[1]=="--resume"
        @test parse(Int,fields[4])==6
        @test parse(Float64,fields[5])==0.93
        @test fields[7]=="true"
        # New accepted iterate is incremental; do not corrupt prior row.
        x[400]=1.1
        stage8_append_accepted_step(trace,"--resume",5,2,x,residual,pre,200)
        @test length(readlines(trace))==3
        @test split(readlines(trace)[3],",")[7]=="false"
    end
end
