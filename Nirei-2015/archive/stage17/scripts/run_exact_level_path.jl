#!/usr/bin/env julia
# Local only: julia --project=. scripts/run_exact_level_path.jl [T] [N] [seed]
# Three-choice 2015 finite-N firms + exact ex-post nonlinear levels.
# Forecasts still use Appendix A: THIS IS NOT EXACT STOCHASTIC RE.
using Nirei2015
using TOML

function main()
    T=length(ARGS)>=1 ? parse(Int,ARGS[1]) : 3
    N=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 48
    seed=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 20261101
    T>0 && N>0 || error("T and N must be positive")
    state=initialize_simulation(;N=N,base_seed=seed)
    records=NamedTuple[]
    failure="none"
    for _ in 1:T
        try
            rec=simulate_exact_level_period!(state)
            push!(records,rec)
            println("t=$(rec.t) K=$(rec.K_next) X=$(rec.X) X_e=$(rec.expected_X) mu=$(rec.mu) stable_gap=$(rec.stable_eigenvector_gap)")
        catch err
            err isa InterruptException && rethrow()
            failure=sprint(showerror,err)
            println(stderr,"EXACT-LEVEL PATH STOP at t=$(state.period): $failure")
            break
        end
    end
    folder=normpath(joinpath(@__DIR__,"..","results"))
    mkpath(folder)
    stem="exact_levels_T$(T)_N$(N)_seed$(seed)"
    csvfile=joinpath(folder,stem*".csv")
    open(csvfile,"w") do io
        headers=isempty(records) ? () : propertynames(first(records))
        println(io,join(string.(headers),","))
        for rec in records
            println(io,join((string(getproperty(rec,k)) for k in headers),","))
        end
    end
    info=Dict(
        "model"=>"2015 original 3-choice finite-N firms, exact ex-post household levels, APPROXIMATE Appendix A forecasts; NOT fully solved stochastic RE",
        "requested_periods"=>T,"completed_periods"=>length(records),
        "completed"=>length(records)==T,"failure"=>failure,
        "N"=>N,"seed"=>seed,
        "max_abs_resource_residual"=>isempty(records) ? 0.0 : maximum(r->abs(r.resources_residual),records),
        "max_abs_labor_residual"=>isempty(records) ? 0.0 : maximum(r->abs(r.labor_residual),records),
        "max_abs_capital_residual"=>isempty(records) ? 0.0 : maximum(r->abs(r.firm_capital_residual),records),
        "max_abs_investment_residual"=>isempty(records) ? 0.0 : maximum(r->abs(r.firm_investment_residual),records),
        "max_abs_stable_gap"=>isempty(records) ? 0.0 : maximum(r->abs(r.stable_eigenvector_gap),records),
    )
    tomlfile=joinpath(folder,stem*".toml")
    open(tomlfile,"w") do io
        TOML.print(io,info)
    end
    println("SUMMARY: ",info)
    println("CSV: ",csvfile," TOML: ",tomlfile)
    length(records)==T || error("exact-level path incomplete at t=$(state.period)")
end
main()
