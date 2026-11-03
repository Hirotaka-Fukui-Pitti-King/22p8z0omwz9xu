# Economic comparison from already converged 200-date Julia market roots.
# No MATLAB fixtures, no GE re-solving, no additional Julia dependencies.
#
# julia --project=. scripts/report_polished_economic_results.jl
#
# CAUTION: The four equilibria are finite-horizon market roots with
# tested Euler/KKT accuracy. Terminal continuation and grid convergence
# are NOT yet certified. The report and data are explicitly PROVISIONAL.

using FMR2020

const ECONOMIC_ROOT_LABELS = (
    boom="credit_expansion_euler_polished",
    nolev="flexible_unleveraged_euler_polished",
    flexible="flexible_leveraged_euler_polished",
    fixed="fixed_leveraged_euler_polished",
)

"""
    economic_read_panel(path)

Read the saved numeric Stage3 market-root CSV with strict row and
column validation. A panel returned from this reader is ONLY a stored
candidate path; reading it does not certify equilibrium.
"""
function economic_read_panel(path::AbstractString)
    isfile(path) || throw(ArgumentError("missing Julia GE root CSV: "*path))
    lines=readlines(path)
    length(lines)>41 ||
        throw(DimensionMismatch("economic transition must exceed the first 40 quarters"))
    names=strip.(split(lines[1],","))
    wanted=("period","B","K","L","Y","I","q","r","w","terms",
            "foreign_consumption","home_consumption","wealth",
            "nominal_exchange_rate")
    length(unique(names))==length(names) ||
        throw(ArgumentError("duplicate CSV headers"))
    all(name->name in names,wanted) ||
        throw(ArgumentError("saved GE root lacks required economic variables"))
    panel=Dict{String,Vector{Float64}}(
        name=>Vector{Float64}(undef,length(lines)-1) for name in wanted)
    for t in 1:length(lines)-1
        fields=strip.(split(lines[t+1],","))
        length(fields)==length(names) ||
            throw(DimensionMismatch("CSV row $t has wrong column count"))
        for name in wanted
            k=findfirst(==(name),names)
            v=parse(Float64,fields[k])
            isfinite(v) || throw(DomainError(v,"nonfinite $name at date $t"))
            panel[name][t]=v
        end
        panel["period"][t]==t ||
            throw(ArgumentError("nonconsecutive economic dates at row $t"))
        for name in ("K","L","Y","q","terms",
                     "foreign_consumption","home_consumption")
            panel[name][t]>0 ||
                throw(DomainError(panel[name][t],
                    "nonpositive $name at date $t"))
        end
    end
    return panel
end

"""
    economic_validate_shared_experiment(boom, nolev, flexible, fixed;
                                        shock_start=41)

Verify every sudden-stop equilibrium was computed for the SAME foreign
credit experiment and inherited physical capital. This does NOT by
itself verify the common *wealth distribution*, which is not stored in
the aggregate CSV and must be established by the generating solver.
"""
function economic_validate_shared_experiment(
    boom, nolev, flexible, fixed; shock_start::Int=41)

    n=length(boom["B"])
    2<=shock_start<=n ||
        throw(ArgumentError("shock date outside boom path"))
    for (name,panel) in (("unleveraged",nolev),
                         ("leveraged flexible",flexible),
                         ("leveraged fixed",fixed))
        length(panel["B"])==n ||
            throw(DimensionMismatch("$name has a different horizon"))
        all(isapprox.(panel["B"],nolev["B"];rtol=1e-10,atol=1e-10)) ||
            throw(ArgumentError("$name has a different foreign-credit shock"))
        isapprox(panel["K"][1],boom["K"][shock_start];
                 rtol=1e-9,atol=1e-9) ||
            throw(ArgumentError("$name has inconsistent inherited K at surprise"))
    end
    return true
end

"""
    economic_log_composite(panel, t, chi)

Log Cobb-Douglas composite consumption, up to an irrelevant shared
normalization constant. All log DIFFERENCES are invariant to this
normalization. Here C_H and C_F are physical home and foreign goods.
"""
function economic_log_composite(panel,t::Int,chi::Real)
    0<chi<1 || throw(ArgumentError("Cobb-Douglas weight must be in (0,1)"))
    return chi*log(panel["home_consumption"][t])+
           (1-chi)*log(panel["foreign_consumption"][t])
end

"""
    economic_response_rows(boom, nolev, flexible, fixed;
                           shock_start=41, horizon=30, chi=0.6)

Two DIFFERENT counterfactual objects:
(1) the *no-surprise continuing boom* at each same calendar date,
and (2) the exchange-rate or leverage contrast holding the credit
reversal and initial conditions fixed.

Response numbers are 100 times LOG ratios (approximately log points).
They are not the paper's welfare compensating variation.
"""
function economic_response_rows(
    boom,nolev,flexible,fixed;
    shock_start::Int=41,horizon::Int=30,chi::Float64=0.6)

    horizon>=1 || throw(ArgumentError("comparison horizon must be positive"))
    length(boom["B"])>=shock_start+horizon-1 ||
        throw(DimensionMismatch("boom no-shock counterfactual too short"))
    all(length(p["B"])>=horizon for p in (nolev,flexible,fixed)) ||
        throw(DimensionMismatch("sudden-stop equilibrium too short"))
    economic_validate_shared_experiment(
        boom,nolev,flexible,fixed;shock_start=shock_start)
    rows=NamedTuple[]
    for h in 0:horizon-1
        t=h+1
        bt=shock_start+h
        cb=economic_log_composite(boom,bt,chi)
        cn=economic_log_composite(nolev,t,chi)
        cf=economic_log_composite(flexible,t,chi)
        cx=economic_log_composite(fixed,t,chi)
        push!(rows,(
            quarter=h,
            boom_calendar_period=bt,
            composite_no_leverage_shock_pct=100*(cn-cb),
            composite_leveraged_flexible_shock_pct=100*(cf-cb),
            composite_leveraged_fixed_shock_pct=100*(cx-cb),
            composite_fixed_minus_flexible_pct=100*(cx-cf),
            composite_leverage_minus_unleveraged_pct=100*(cf-cn),
            home_fixed_minus_flexible_pct=100*
                log(fixed["home_consumption"][t]/
                    flexible["home_consumption"][t]),
            foreign_fixed_minus_flexible_pct=100*
                log(fixed["foreign_consumption"][t]/
                    flexible["foreign_consumption"][t]),
            output_fixed_minus_flexible_pct=100*
                log(fixed["Y"][t]/flexible["Y"][t]),
            investment_fixed_minus_flexible_level=
                fixed["I"][t]-flexible["I"][t],
            terms_fixed_minus_flexible_pct=100*
                log(fixed["terms"][t]/flexible["terms"][t]),
            wage_fixed_minus_flexible_pct=100*
                log(fixed["w"][t]/flexible["w"][t]),
            rate_fixed_minus_flexible_pp=100*
                (fixed["r"][t]-flexible["r"][t]),
        ))
    end
    return rows
end

function economic_write_rows(rows,destination::AbstractString)
    isempty(rows) && throw(ArgumentError("no economic comparison rows"))
    keys0=propertynames(first(rows))
    mkpath(dirname(destination))
    open(destination,"w") do io
        println(io,join(string.(keys0),","))
        for row in rows
            propertynames(row)==keys0 ||
                throw(ArgumentError("inconsistent comparison columns"))
            println(io,join((string(getproperty(row,k)) for k in keys0),","))
        end
    end
    return destination
end

function economic_summary(rows,boom,flexible,fixed,chi;
                          shock_start::Int=41)
    n=length(rows)
    n>0 || throw(ArgumentError("empty comparison"))
    start=first(rows)
    cumulative_c_flex=0.0
    cumulative_c_fixed=0.0
    cumulative_c_boom=0.0
    for h in 0:n-1
        t=h+1
        bt=shock_start+h
        cumulative_c_flex+=exp(economic_log_composite(flexible,t,chi))
        cumulative_c_fixed+=exp(economic_log_composite(fixed,t,chi))
        cumulative_c_boom+=exp(economic_log_composite(boom,bt,chi))
    end
    average_fx= sum(r.composite_fixed_minus_flexible_pct for r in rows)/n
    max_fx=maximum(abs(r.composite_fixed_minus_flexible_pct) for r in rows)
    max_lev=maximum(abs(r.composite_leverage_minus_unleveraged_pct) for r in rows)
    return (
        impact_flexible_vs_no_shock=start.composite_leveraged_flexible_shock_pct,
        impact_fixed_vs_no_shock=start.composite_leveraged_fixed_shock_pct,
        impact_fixed_minus_flexible=start.composite_fixed_minus_flexible_pct,
        impact_leverage_minus_unleveraged=start.composite_leverage_minus_unleveraged_pct,
        average_fixed_minus_flexible=average_fx,
        maxabs_fixed_minus_flexible=max_fx,
        maxabs_leverage_minus_unleveraged=max_lev,
        cumulative_flexible_vs_boom_pct=100*
            (cumulative_c_flex/cumulative_c_boom-1),
        cumulative_fixed_vs_boom_pct=100*
            (cumulative_c_fixed/cumulative_c_boom-1),
        cumulative_fixed_minus_flexible_pct=100*
            (cumulative_c_fixed/cumulative_c_flex-1),
    )
end

function economic_write_summary(summary,destination,horizon::Int)
    mkpath(dirname(destination))
    open(destination,"w") do io
        println(io,"# PROVISIONAL economic comparison — Euler-polished Stage3")
        println(io)
        println(io,"**Not certified as infinite-horizon GE.** ",
                "Household Euler/KKT and 599-equation market residuals ",
                "converged on the stored T=200 paths. Terminal budget ",
                "and grid/horizon convergence for all shock regimes ",
                "remain unresolved.")
        println(io)
        println(io,"Comparison window: $horizon quarters from the ",
                "unexpected shock, with the time-matched continuing boom ",
                "as the no-surprise counterfactual.")
        println(io)
        println(io,"Log-response values below are 100 × differences in ",
                "log composite consumption. Positive fixed-minus-flexible ",
                "means consumption is HIGHER under fixed exchange rates.")
        println(io)
        println(io,"| Economic contrast | Value |")
        println(io,"|---|---:|")
        for (label,key) in (
            ("On-impact leveraged flexible vs no surprise",
             :impact_flexible_vs_no_shock),
            ("On-impact leveraged fixed vs no surprise",
             :impact_fixed_vs_no_shock),
            ("On-impact fixed minus flexible",
             :impact_fixed_minus_flexible),
            ("On-impact leverage effect under flexible",
             :impact_leverage_minus_unleveraged),
            ("Average fixed minus flexible ($horizon quarters)",
             :average_fixed_minus_flexible),
            ("Largest absolute fixed-versus-flexible log gap",
             :maxabs_fixed_minus_flexible),
            ("Largest absolute leverage log gap",
             :maxabs_leverage_minus_unleveraged),
            ("Cumulative consumption level flexible vs continuing boom",
             :cumulative_flexible_vs_boom_pct),
            ("Cumulative consumption level fixed vs continuing boom",
             :cumulative_fixed_vs_boom_pct),
            ("Cumulative consumption level fixed minus flexible",
             :cumulative_fixed_minus_flexible_pct),
        )
            println(io,"| ",label," | ",
                    round(getproperty(summary,key);digits=7)," |")
        end
        println(io)
        println(io,"**Interpretation limits.** The cumulative entries ",
                "compare unweighted quarter sums of Cobb-Douglas ",
                "composite consumption indexes; they are NOT household ",
                "welfare or compensating variation. No assumption ",
                "about terminal welfare is imposed in this report.")
    end
    return destination
end

function report_polished_main(root::AbstractString=
    normpath(joinpath(@__DIR__,"..","results","economic_reproduction")))

    files=(joinpath(root,"economic_"*ECONOMIC_ROOT_LABELS.boom*
                     "_MARKET_ROOT_UNCERTIFIED.csv"),
           joinpath(root,"economic_"*ECONOMIC_ROOT_LABELS.nolev*
                     "_MARKET_ROOT_UNCERTIFIED.csv"),
           joinpath(root,"economic_"*ECONOMIC_ROOT_LABELS.flexible*
                     "_MARKET_ROOT_UNCERTIFIED.csv"),
           joinpath(root,"economic_"*ECONOMIC_ROOT_LABELS.fixed*
                     "_MARKET_ROOT_UNCERTIFIED.csv"))
    boom,nolev,flexible,fixed=map(economic_read_panel,files)
    chi=Float64(FMR2020.default_parameters().chi)
    rows=economic_response_rows(boom,nolev,flexible,fixed;
                                shock_start=41,horizon=30,chi=chi)
    s=economic_summary(rows,boom,flexible,fixed,chi;shock_start=41)
    outputrows=joinpath(root,"PROVISIONAL_figure4_5_economic_comparison.csv")
    outputbrief=joinpath(root,"PROVISIONAL_economic_comparison.md")
    economic_write_rows(rows,outputrows)
    economic_write_summary(s,outputbrief,30)
    println("PROVISIONAL ECONOMIC RESULTS: ",s)
    println("Quarter-by-quarter economic comparison: ",outputrows)
    println("Research-summary table: ",outputbrief)
    println("NOT a certified infinite-horizon welfare comparison.")
    return s
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)<=1 ||
        error("Usage: julia --project=. scripts/report_polished_economic_results.jl [RESULTS_DIR]")
    root=isempty(ARGS) ?
         normpath(joinpath(@__DIR__,"..","results","economic_reproduction")) :
         ARGS[1]
    report_polished_main(root)
end
