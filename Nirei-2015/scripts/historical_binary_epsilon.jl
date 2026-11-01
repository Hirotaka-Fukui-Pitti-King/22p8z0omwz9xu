#!/usr/bin/env julia
# Historically source-grounded 2013 combined protocol: one-sided binary
# indivisible investment AND continuum capital Eq (18) epsilon iteration.
# For 2015 comparison: published Table 1 sigma=1.5, no altered shocks.
# This is NOT the published 2015 three-choice finite-N CES fixed point;
# we explicitly reaggregate firms and measure its exact CES residual.
# julia --project=. scripts/historical_binary_epsilon.jl [T] [N] [seed]

using Nirei2015
using TOML
using Statistics: std, cor

function historical_binary_period!(state::SimulationState)
    cal = state.calibration
    p = cal.parameters
    s = cal.steady_state
    firm = state.firms
    n = length(firm.capital)
    invn = inv(Float64(n))
    prior_K = state.K
    prior_mu = state.mu_prior
    expected = compute_expected_equilibrium(state.re,cal,prior_K,prior_mu)
    draw_productivity!(firm,state.rngs,cal)
    d = Vector{Float64}(undef,n)
    threshold = Vector{Float64}(undef,n)
    gain = Vector{Float64}(undef,n)
    @inbounds for i in 1:n
        d[i]=(1.0-p.delta)*firm.capital[i]
        lambda=firm.lambda[i]
        b=threshold_coefficient(lambda,firm.next_productivity[i],cal)
        threshold[i]=d[i]/b
        gain[i]=d[i]*(lambda-1.0)*invn
    end
    order=sortperm(threshold)
    baseline=(1.0-p.delta)*prior_K
    slope=cal.threshold.A^(cal.constants.rho-1.0)/cal.constants.rho
    # Before the first switch no firms invest; X starts at zero.
    total_X=0.0
    left=0.0
    index=1
    selected_K=NaN
    selected_X=NaN
    selected_score=Inf
    n_roots=0
    n_positive=0
    while index<=n
        upper=threshold[order[index]]
        Khat=baseline+slope*total_X
        if left<Khat<=upper
            n_roots+=1
            if total_X>0.0
                n_positive+=1
                score=abs(log(total_X/expected.X_next))
                if score<selected_score
                    selected_K=Khat
                    selected_X=total_X
                    selected_score=score
                end
            end
        end
        while index<=n && threshold[order[index]]==upper
            i=order[index]
            total_X+=gain[i]
            index+=1
        end
        left=upper
    end
    Khat=baseline+slope*total_X
    if Khat>left
        n_roots+=1
        if total_X>0.0
            n_positive+=1
            score=abs(log(total_X/expected.X_next))
            if score<selected_score
                selected_K=Khat
                selected_X=total_X
                selected_score=score
            end
        end
    end
    isfinite(selected_K) && selected_X>0.0 ||
        error("No positive-X historical binary continuum equilibrium at t=$(state.period): roots=$n_roots positive=$n_positive")
    next_k=Vector{Float64}(undef,n)
    direct_X=0.0
    investing=0
    @inbounds for i in 1:n
        is_investing=selected_K>threshold[i]
        next_k[i]=is_investing ? d[i]*firm.lambda[i] : d[i]
        direct_X+=next_k[i]-d[i]
        investing+=is_investing
    end
    direct_X*=invn
    abs(direct_X-selected_X)<1e-9*max(1.0,selected_X) ||
        error("historical one-sided individual investment mismatch: $(direct_X-selected_X)")
    true_CES_K=aggregate_capital(firm.next_productivity,next_k,p.alpha,cal.constants.rho)
    CES_gap=true_CES_K-selected_K
    firm.capital=next_k
    advance_productivity!(firm)
    logKprior=log(prior_K/s.K)
    logmu=state.re.h*log(true_CES_K/s.K)
    next_mu=s.mu*exp(logmu)
    m=Nirei2015._realized_macroeconomy(
        cal,state.re,logKprior,log(selected_X/s.X),logmu,log(prior_mu/s.mu))
    resource_gap=(m.Y-m.C-selected_X)/m.Y
    state.K=true_CES_K
    state.mu_prior=next_mu
    state.period+=1
    return (t=state.period-1,K_CES=true_CES_K,K_continuum=selected_K,
            CES_gap=CES_gap,X=selected_X,X_expected=expected.X_next,
            epsilon=log(selected_X/expected.X_next),n_roots=n_roots,
            n_positive=n_positive,n_investing=investing,
            direct_X_gap=direct_X-selected_X,Y=m.Y,C=m.C,L=m.L,
            resource_gap=resource_gap)
end

function main()
    T=length(ARGS)>=1 ? parse(Int,ARGS[1]) : 20
    N=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 350000
    seed=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 20261101
    T>=1 && N>=1 || error("T and N must be positive")
    state=initialize_simulation(N=N,base_seed=seed)
    records=NamedTuple[]
    failure="none"
    for _ in 1:T
        try
            rec=historical_binary_period!(state)
            push!(records,rec)
            if rec.t<=5 || rec.t==T-1
                println("t=$(rec.t) K_continuum=$(rec.K_continuum) K_CES=$(rec.K_CES) CES_gap=$(rec.CES_gap) X=$(rec.X) Xe=$(rec.X_expected) epsilon=$(rec.epsilon) roots=$(rec.n_roots)")
            end
        catch err
            failure=sprint(showerror,err)
            println(stderr,"HISTORICAL BINARY+EPSILON STOP at t=$(state.period): ",failure)
            break
        end
    end
    dest=joinpath(@__DIR__,"..","results")
    mkpath(dest)
    prefix="historical_binary_epsilon_T$(T)_N$(N)_seed$(seed)"
    csv=joinpath(dest,prefix*".csv")
    cols=(:t,:K_CES,:K_continuum,:CES_gap,:X,:X_expected,:epsilon,:n_roots,
          :n_positive,:n_investing,:direct_X_gap,:Y,:C,:L,:resource_gap)
    open(csv,"w") do io
        println(io,join(string.(cols),","))
        for r in records
            println(io,join((string(getfield(r,k)) for k in cols),","))
        end
    end
    info=Dict(
        "description"=>"Historical 2013 one-sided binary + continuum Eq18 epsilon protocol. NOT published 2015 three-choice exact CES GE.",
        "T"=>T,"N"=>N,"seed"=>seed,"completed"=>length(records)==T,
        "completed_periods"=>length(records),"failure"=>failure,
        "max_abs_CES_gap"=>isempty(records) ? 0.0 : maximum(r->abs(r.CES_gap),records),
        "max_abs_X_gap"=>isempty(records) ? 0.0 : maximum(r->abs(r.direct_X_gap),records),
        "max_abs_resource_gap"=>isempty(records) ? 0.0 : maximum(r->abs(r.resource_gap),records))
    if length(records)>100
        r=records[101:end]
        info["postburnin"]=Dict(
            "Y_sd_percent"=>100*std(log(z.Y) for z in r),
            "C_sd_percent"=>100*std(log(z.C) for z in r),
            "X_sd_percent"=>100*std(log(z.X) for z in r),
            "L_sd_percent"=>100*std(log(z.L) for z in r),
            "K_sd_percent"=>100*std(log(z.K_CES) for z in r),
            "Y_X_correlation"=>cor([log(z.Y) for z in r],[log(z.X) for z in r]))
    end
    output=joinpath(dest,prefix*".toml")
    open(output,"w") do io
        TOML.print(io,info)
    end
    println("RESULT: ",info)
    println("Saved: ",csv," / ",output)
    length(records)==T || error("historical binary epsilon incomplete at t=$(state.period)")
end

main()
