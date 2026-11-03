using Test

include(joinpath(@__DIR__,"..","scripts","report_polished_economic_results.jl"))

function synthetic_economic_panel(n::Int;
                                  home::Float64=2.0,
                                  foreign::Float64=1.0)
    n>=80 || throw(ArgumentError("test requires a 30-quarter surprise window"))
    p=Dict{String,Vector{Float64}}()
    names=("period","B","K","L","Y","I","q","r","w","terms",
           "foreign_consumption","home_consumption","wealth",
           "nominal_exchange_rate")
    for name in names
        p[name]=fill(1.0,n)
    end
    p["period"]=Float64.(1:n)
    p["B"]=fill(-4.0,n)
    p["K"]=fill(8.0,n)
    p["I"]=fill(0.2,n)
    p["home_consumption"]=fill(home,n)
    p["foreign_consumption"]=fill(foreign,n)
    return p
end

function synthetic_economic_csv(path::String,p)
    header=("period","B","K","L","Y","I","q","r","w","terms",
            "foreign_consumption","home_consumption","wealth",
            "nominal_exchange_rate")
    open(path,"w") do io
        println(io,join(header,","))
        for t in eachindex(p["period"])
            println(io,join((string(p[name][t]) for name in header),","))
        end
    end
end

@testset "Economic fixed-versus-flexible contrast from saved market roots" begin
    boom=synthetic_economic_panel(80)
    nolev=synthetic_economic_panel(80;home=1.7,foreign=0.85)
    flexible=synthetic_economic_panel(80;home=1.6,foreign=0.8)
    fixed=synthetic_economic_panel(80;home=1.8,foreign=0.9)
    @test economic_validate_shared_experiment(boom,nolev,flexible,fixed)

    rows=economic_response_rows(boom,nolev,flexible,fixed;
                                shock_start=41,horizon=30,chi=0.6)
    @test length(rows)==30
    @test rows[1].quarter==0
    @test rows[end].boom_calendar_period==70
    @test isapprox(rows[1].composite_leveraged_flexible_shock_pct,
                   100log(0.8);atol=1e-12)
    @test isapprox(rows[1].composite_leveraged_fixed_shock_pct,
                   100log(0.9);atol=1e-12)
    @test isapprox(rows[1].composite_fixed_minus_flexible_pct,
                   100log(0.9/0.8);atol=1e-12)
    @test isapprox(rows[1].composite_leverage_minus_unleveraged_pct,
                   100log(1.6/1.7);atol=1e-12)

    result=economic_summary(rows,boom,flexible,fixed,0.6)
    @test result.impact_fixed_minus_flexible>0
    @test isapprox(result.cumulative_flexible_vs_boom_pct,-20.0;atol=1e-10)
    @test isapprox(result.cumulative_fixed_vs_boom_pct,-10.0;atol=1e-10)
    @test isapprox(result.cumulative_fixed_minus_flexible_pct,12.5;atol=1e-10)

    bad=deepcopy(fixed)
    bad["B"][10]+=0.1
    @test_throws ArgumentError economic_response_rows(
        boom,nolev,flexible,bad;chi=0.6)
    bad=deepcopy(fixed)
    bad["K"][1]+=0.1
    @test_throws ArgumentError economic_validate_shared_experiment(
        boom,nolev,flexible,bad)

    mktempdir() do dir
        paths=(joinpath(dir,"economic_"*ECONOMIC_ROOT_LABELS.boom*
                         "_MARKET_ROOT_UNCERTIFIED.csv"),
               joinpath(dir,"economic_"*ECONOMIC_ROOT_LABELS.nolev*
                         "_MARKET_ROOT_UNCERTIFIED.csv"),
               joinpath(dir,"economic_"*ECONOMIC_ROOT_LABELS.flexible*
                         "_MARKET_ROOT_UNCERTIFIED.csv"),
               joinpath(dir,"economic_"*ECONOMIC_ROOT_LABELS.fixed*
                         "_MARKET_ROOT_UNCERTIFIED.csv"))
        for (path,p) in zip(paths,(boom,nolev,flexible,fixed))
            synthetic_economic_csv(path,p)
            parsed=economic_read_panel(path)
            @test parsed["B"]==p["B"]
            @test parsed["home_consumption"]==p["home_consumption"]
        end
        summary=report_polished_main(dir)
        @test isapprox(summary.impact_fixed_minus_flexible,
                       100log(0.9/0.8);atol=1e-12)
        csvout=joinpath(dir,"PROVISIONAL_figure4_5_economic_comparison.csv")
        mdout=joinpath(dir,"PROVISIONAL_economic_comparison.md")
        @test isfile(csvout)
        @test length(readlines(csvout))==31
        @test isfile(mdout)
        @test occursin("Not certified",read(mdout,String))
        @test_throws ArgumentError economic_read_panel(
            joinpath(dir,"unknown.csv"))
    end
end
