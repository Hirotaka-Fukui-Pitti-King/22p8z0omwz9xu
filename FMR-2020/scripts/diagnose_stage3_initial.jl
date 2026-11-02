# Stage 3 initial-equilibrium diagnostic. Deliberately evaluates each 399-equation
# residual ONCE and does not launch nonlinear root-finding.
# See original Main.m lines 127-154 and 342-351: no-leverage first.
using FMR2020
using LinearAlgebra

function stage3_initial_guess(steady, context)
    inp=context.economy
    n=length(inp.B)
    x=Vector{Float64}(undef,2*n-1)
    for t in 1:n-1
        x[t]=steady.calibrated.K+
             (inp.initial.K-steady.calibrated.K)*0.8^t
    end
    x[n:end].=steady.calibrated.L
    return x
end

function stage3_boom_guess(filename, inp)
    n=length(inp.B)
    if !isfile(filename)
        return FMR2020.initial_transition_guess(inp)
    end
    lines=readlines(filename)
    length(lines)==n+1 ||
        error("Stage 2 warm-start CSV has unexpected horizon")
    head=split(lines[1],",")
    k=findfirst(==("K"),head)
    l=findfirst(==("L"),head)
    (isnothing(k) || isnothing(l)) &&
        error("Stage 2 warm-start CSV missing K/L columns")
    K=Float64[]
    L=Float64[]
    for row in lines[2:end]
        fields=split(row,",")
        push!(K,parse(Float64,fields[k]))
        push!(L,parse(Float64,fields[l]))
    end
    return vcat(K[2:end],L)
end

function main()
    root=normpath(joinpath(@__DIR__,".."))
    steady=solve_three_steady_states(default_parameters();verbose=false)
    boomctx=stage3_expansion_context(steady)
    warm=stage3_boom_guess(joinpath(root,"results",
                                    "stage2_transition_expansion.csv"),
                           boomctx.economy)
    boom=solve_stage3_flexible(boomctx;x0=warm,maxiter=3,
                               tolerance=1e-6,verbose=false)
    println("STAGE3 INITIAL AUDIT: corrected foreign boom residual=",
            boom.maximum_scaled_residual)
    for leverage in (0.0,1.0/16.0)
        ctx=stage3_reversal_context(steady,boom;
                                    shock_start=41,
                                    leverage_hat=leverage)
        x=stage3_initial_guess(steady,ctx)
        work=FMR2020.TransitionWorkspace(ctx.economy)
        residual=zeros(Float64,length(x))
        println("=== SURPRISE REVERSAL leverage=",leverage," ===")
        try
            FMR2020.stage3_residual!(residual,x,ctx,work)
            n=length(ctx.economy.B)
            aview=view(residual,1:n-1)
            lview=view(residual,n:2*n-1)
            println("total norm=",norm(residual),
                    " max asset=",maximum(abs, aview),
                    " max labor=",maximum(abs,lview))
            println("labor residual periods 1-5=",collect(lview[1:5]))
            println("asset residual periods 1-5=",collect(aview[1:5]))
            println("initial terms of trade=",work.terms[1],
                    " pre-shock terms=",ctx.terms_pre)
            println("initial C_H=",work.aggregate_home[1],
                    " wage=",work.wage[1],
                    " hours=",work.L[1],
                    " implied hours=",work.labor_implied[1])
            println("mass check=",sum(view(work.distribution,:,:,1)),
                    " FX initial gross return=",1+work.rate[1])
        catch err
            println("INITIAL RESIDUAL INFEASIBLE: ",sprint(showerror,err))
            rethrow()
        end
    end
    println("DIAGNOSTIC COMPLETE: no full equilibrium solved in this script")
end

main()
