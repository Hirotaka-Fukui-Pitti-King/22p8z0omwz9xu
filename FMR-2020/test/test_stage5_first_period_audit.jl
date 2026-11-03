using Test
using FMR2020

include(joinpath(@__DIR__,"..","scripts","stage5_initial_fx_audit.jl"))

@testset "Stage 5 initial FX household portfolio accounting and branches" begin
    rows=stage5_audit_main()
    @test length(rows)==3
    for row in rows
        @test row.maximum_scaled_residual<1e-6
        @test isapprox(row.qK+row.B_initial,
                       row.pre_distribution_wealth;atol=2e-5)
        @test isapprox(row.predicted_consumption_gap,
                       row.actual_consumption_gap;atol=5e-5,rtol=1e-2)
        @test row.lower_trade_root.root<row.trade_price_peak
        @test row.upper_trade_root!==nothing
        @test row.upper_trade_root.root>row.trade_price_peak
        @test isfinite(row.observed_aggregate_scalar_root_gap)
    end
    # A conditional trade-price root is NOT an independent GE root.
    # The observed FX path lies close to the higher branch; the
    # lower root implies a depreciation but can violate budgets.
    @test rows[2].lower_trade_root.root<rows[2].terms_pre
    @test rows[2].upper_trade_root.root>rows[2].terms_pre
    @test rows[2].terms_observed>rows[2].terms_pre
    @test abs(rows[1].gross_capital_gap)<1e-4
    @test abs(rows[2].gross_capital_gap)>1e-3
    @test abs(rows[3].gross_capital_gap)>1e-3
end
