using Test
using FMR2020

include(joinpath(@__DIR__, "..", "scripts", "stage8_compare_saved_roots.jl"))

@testset "Stage 8 independent terminal-budget algebra" begin
    # Actual and terminal steady-state budget at the SAME wealth
    # and stationary continuation (c,a'), with independently varied prices.
    a=0.65
    productivity=1.1
    c=0.8
    saving=0.2
    r=0.015
    wl=1.7
    profits=0.14
    m=1.65
    rs=0.01
    wls=1.68
    profitss=0.13
    ms=1.64
    z=stage8_terminal_budget_components(a,productivity,c,saving,
                                          r,wl,profits,m,rs,wls,profitss,ms)
    @test abs(z.decomposition_error)<1e-14
    @test isapprox(z.interest_contribution,(r-rs)*a;atol=1e-14)
    @test isapprox(z.wage_contribution,(wl-wls)*productivity;atol=1e-14)
    @test isapprox(z.profit_contribution,profits-profitss;atol=1e-14)
    @test isapprox(z.home_cost_contribution,-(m-ms)*c;atol=1e-14)
    @test isapprox(
        z.actual_budget,
        z.stationary_policy_budget+z.interest_contribution+
        z.wage_contribution+z.profit_contribution+z.home_cost_contribution;
        atol=1e-14)
end

@testset "Stage 8 identical-price continuation matches stationary budget" begin
    r=0.01
    wl=1.2
    profits=.2
    m=1.7
    a=1.5
    z=1.2
    save=.35
    c=((1+r)*a+wl*z+profits-save)/m
    parts=stage8_terminal_budget_components(a,z,c,save,
                                              r,wl,profits,m,
                                              r,wl,profits,m)
    @test abs(parts.actual_budget)<1e-14
    @test abs(parts.stationary_policy_budget)<1e-14
    @test parts.interest_contribution==0
    @test parts.wage_contribution==0
    @test parts.profit_contribution==0
    @test parts.home_cost_contribution==0
    @test abs(parts.decomposition_error)<1e-14
end

@testset "Stage 8 rejected nonfinite comparison data" begin
    @test_throws ArgumentError stage8_terminal_budget_components(
        NaN,1.0,1.0,0.1,0.01,1.0,0.1,1.0,0.01,1.0,0.1,1.0)
end

@testset "Stage 8 terminal budget probability-mass diagnostics without clipping" begin
    tiny=nextfloat(0.0)
    masses=[1.0-1e-6,1e-6,tiny,0.0]
    budgets=[2e-4,-0.005,-0.09,NaN]
    results=stage8_terminal_mass_summary(masses,budgets)
    @test results.material_cell_count==2
    @test results.material_cutoff==1e-12
    @test results.max_budget_all_occupied==0.09
    @test results.max_budget_material_cells==0.005
    @test isapprox(results.material_probability_mass,1.0;atol=1e-15)
    @test isapprox(results.violation_mass_above_1e3,1e-6;atol=1e-15)
    @test isapprox(results.weighted_absolute_budget,
                   (1.0-1e-6)*2e-4+1e-6*0.005;atol=1e-15)
    @test_throws ArgumentError stage8_terminal_mass_summary(
        [-0.1,1.1],[0.1,0.2])
    @test_throws ArgumentError stage8_terminal_mass_summary(
        [0.5,0.5],[0.1,NaN])
    @test_throws DimensionMismatch stage8_terminal_mass_summary(
        [1.0],[0.2,0.3])
end

println("STAGE8 TWO ROOT COMPARISON LIGHT TESTS PASS; no full GE solver")
