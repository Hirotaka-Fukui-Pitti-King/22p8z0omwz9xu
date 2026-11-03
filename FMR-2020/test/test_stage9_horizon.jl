using Test
using FMR2020

include(joinpath(@__DIR__,"..","scripts","stage9_horizon_extension.jl"))

@testset "Stage 9 strict exact shock/horizon invariants" begin
    paths=stage9_filepaths(normpath(joinpath(@__DIR__,"..")))
    steady=FMR2020.solve_three_steady_states(
        FMR2020.default_parameters();verbose=false)
    boom=stage5_existing_boom(steady,paths.prior)
    old,long,inv=stage9_contexts(steady,boom)
    @test inv.old_horizon==200
    @test inv.new_horizon==300
    @test inv.old_unknown_count==599
    @test inv.new_unknown_count==899
    @test inv.same_terminal
    @test inv.same_initial_distribution
    @test long.economy.B[1:200]==old.economy.B
    @test all(==(old.economy.terminal.B),long.economy.B[201:300])
    @test old.mu_pre==long.mu_pre
    @test isequal(old.economy.initial,long.economy.initial)
    @test isequal(old.economy.terminal,long.economy.terminal)
    @test isapprox(sum(long.mu_pre),1.0;atol=1e-10)
    @test maximum(abs,long.economy.B[1:200]-old.economy.B)==0
    @test maximum(abs,long.economy.z[1:200]-old.economy.z)==0
    @test_throws DimensionMismatch FMR2020.stage3_prolong_aggregate_guess(
        zeros(598),old,long)
    @test FMR2020.stage3_extend_horizon(old,200)===old
    @test_throws ArgumentError FMR2020.stage3_extend_horizon(old,199)
end

@testset "Stage 9 prolong BOTH saved strict roots without touching old dates" begin
    base=normpath(joinpath(@__DIR__,".."))
    paths=stage9_filepaths(base)
    steady=FMR2020.solve_three_steady_states(
        FMR2020.default_parameters();verbose=false)
    boom=stage5_existing_boom(steady,paths.prior)
    old,long,_=stage9_contexts(steady,boom)
    for (branch,path) in (
        ("appreciation",paths.oldroot.appreciation),
        ("depreciation",paths.oldroot.depreciation),
    )
        original=read_economic_root_guess(path,old)
        @test original!==nothing
        ans=stage9_extend_saved_root(original,old,long)
        @test length(ans)==899
        @test ans[1:199]==original[1:199]
        @test ans[300:499]==original[200:399]
        @test ans[600:799]==original[400:599]
        @test all(isfinite,ans)
        @test all(>(0),ans)
        if branch=="depreciation"
            @test ans[600]<long.terms_pre
        else
            @test ans[600]>long.terms_pre
        end
        @test_throws DimensionMismatch stage9_extend_saved_root(
            original[1:end-1],old,long)
    end
end

println("STAGE9 599-to-899 horizon structural regression tests PASS (no 899 GE root solve)")
