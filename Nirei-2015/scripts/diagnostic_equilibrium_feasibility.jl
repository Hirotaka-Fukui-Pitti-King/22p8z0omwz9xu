#!/usr/bin/env julia
# Nirei (2015) Issue #5, Stage 10. Diagnostic only, NO shock rejection.
# Run from Nirei-2015:
# julia --project=. scripts/diagnostic_equilibrium_feasibility.jl [N] [seed] [draws]
#
# Three independently reconstructed original benchmark dates: 0, 4 and 5.
# Original t=5 productivity profile is analyzed, NEVER replaced. Separate
# Xoshiro draws at its FROZEN state estimate positive-X feasibility.
# This is NOT a nonlinear rational-expectations solution.

using Nirei2015
using Random: Xoshiro, randn!
using Statistics: mean
using TOML

const STAGE10_MC_TAG = 0x2b7e151628aed2a6

function stage10_benchmark_state(t::Int, n::Int, seed::Int)
    state=initialize_simulation(;N=n,base_seed=seed)
    for j in 0:(t-1)
        rec=simulate_period!(state)
        rec.t==j || error("original benchmark timeline changed")
    end
    state.period==t || error("incorrect benchmark date")
    return state
end

# Independent direct firm-by-firm check. Rebuild Eq. (16) coefficient from
# primitives rather than using cached production switches or CES sums.
# The 3-choice policy at a candidate root is kept EXACTLY as specified.
function stage10_direct_response(state, candidate::Float64)
    cal=state.calibration
    p=cal.parameters
    rho=cal.constants.rho
    firms=state.firms
    n=length(firms.capital)
    total_power=0.0
    total_x=0.0
    investing=0
    divesting=0
    inactive=0
    min_threshold_distance=Inf
    @inbounds for i in 1:n
        li=firms.lambda[i]
        ai=firms.next_productivity[i]
        d=(1.0-p.delta)*firms.capital[i]
        z=log(li)
        b=exp(log(cal.threshold.A)+
              ((rho/p.alpha)*log(ai)+log(expm1(rho*z))-log(expm1(z)))/
              (1.0-rho))
        ell=d/(li*b)
        upper=d/b
        policy=candidate<=ell ? -1 : candidate<=upper ? 0 : 1
        investing+=(policy==1)
        divesting+=(policy==-1)
        inactive+=(policy==0)
        mult=li^policy
        knext=d*mult
        total_x+=knext-d
        total_power+=ai^(rho/p.alpha)*knext^rho
        min_threshold_distance=min(
            min_threshold_distance,abs(candidate-ell),abs(candidate-upper)
        )
    end
    K=(total_power/n)^(inv(rho))
    X=total_x/n
    return (
        K=K,X=X,capital_power=total_power/n,
        n_investing=investing,n_divesting=divesting,n_inactive=inactive,
        nearest_threshold_abs_distance=min_threshold_distance
    )
end


# X(K) is monotone in K. At a breakpoint b the old action still holds,
# so positivity starts strictly ABOVE the first crossing to X>0.
# This identifies the exact feasibility cutoff from the original 2N events.
function stage10_positive_investment_cutoff(state)
    ws=state.workspace
    firms=state.firms
    n=length(firms.lambda)
    x=0.0
    @inbounds for i in 1:n
        x+=ws.depreciated[i]*(inv(firms.lambda[i])-1.0)/n
    end
    j=1
    while j<=2n
        boundary=ws.switches[ws.order[j]]
        while j<=2n && ws.switches[ws.order[j]]==boundary
            event=ws.order[j]
            i=(event+1)>>>1
            d=ws.depreciated[i]
            li=firms.lambda[i]
            x+=isodd(event) ?
                d*(1.0-inv(li))/n :
                d*(li-1.0)/n
            j+=1
        end
        if x>0
            return boundary
        end
    end
    error("finite-firm high-capital positive investment was not found")
end

function stage10_wilson(k::Int,n::Int)
    n>0 || error("MC draws must be positive")
    z=1.959963984540054
    phat=k/n
    den=1+z^2/n
    mid=(phat+z^2/(2n))/den
    width=z*sqrt(phat*(1-phat)/n+z^2/(4n^2))/den
    return (low=max(0.0,mid-width),high=min(1.0,mid+width))
end

function stage10_geometry!(state, t::Int;
                           independently_check_every_root::Bool=false)
    cal=state.calibration
    rho=cal.constants.rho
    p=cal.parameters
    priorK=state.K
    priorMu=state.mu_prior
    expected=compute_expected_equilibrium(state.re,cal,priorK,priorMu)
    expected_K18=(1-p.delta)*priorK+
        cal.threshold.A^(rho-1)/rho*expected.X_next
    draw_productivity!(state.firms,state.rngs,cal)
    enumeration=enumerate_fixed_points!(state.workspace,state.firms,cal)
    ws=state.workspace
    count=enumeration.total
    count>0 || error("no capital fixed points at t=$t")
    cutoff=stage10_positive_investment_cutoff(state)
    fixed=stage10_direct_response(state,expected_K18)
    # Note: fixed is a POLICY RESPONSE, not a self-consistent equilibrium.
    nearest=argmin([abs(ws.fixed_capital[j]-expected_K18) for j in 1:count])
    positive=[j for j in 1:count if ws.fixed_investment[j]>0]
    closest_positive=isempty(positive) ? 0 :
        positive[argmin([abs(ws.fixed_capital[j]-expected_K18) for j in positive])]
    selected=isempty(positive) ? 0 :
        select_es1(ws,expected.X_next,expected_K18).index
    ((ws.fixed_capital[count]>cutoff)==!isempty(positive)) ||
        error("break-even capital cutoff contradicts fixed point positivity")
    keys=unique([1,count,nearest,selected,closest_positive])
    check_set=Set(independently_check_every_root ? collect(1:count) :
                  [j for j in keys if j>0])
    rows=Vector{Dict{String,Any}}()
    max_verified_residual=0.0
    max_verified_x_residual=0.0
    for j in 1:count
        Kj=ws.fixed_capital[j]
        Xj=ws.fixed_investment[j]
        row=Dict{String,Any}(
            "t"=>t,"rank"=>j,"K"=>Kj,"X"=>Xj,
            "positive_X"=>Xj>0.0,
            "X_minus_forecast"=>Xj-expected.X_next,
            "logK_minus_forecast"=>log(Kj/expected_K18),
            "nearest_to_expected"=>j==nearest,
            "nearest_positive_to_expected"=>j==closest_positive,
            "ES1_selected"=>j==selected
        )
        if j in check_set
            direct=stage10_direct_response(state,Kj)
            kr=direct.K-Kj
            xr=direct.X-Xj
            max_verified_residual=max(max_verified_residual,abs(kr))
            max_verified_x_residual=max(max_verified_x_residual,abs(xr))
            isapprox(direct.K,Kj;atol=1e-9,rtol=1e-10) ||
                error("direct CES fixed-point root verification failed at t=$t rank=$j")
            isapprox(direct.X,Xj;atol=1e-9,rtol=1e-10) ||
                error("direct investment root verification failed at t=$t rank=$j")
            row["direct_K_residual"]=kr
            row["direct_X_residual"]=xr
            row["investment_share"]=direct.n_investing/length(state.firms.capital)
            row["divestment_share"]=direct.n_divesting/length(state.firms.capital)
            row["inactive_share"]=direct.n_inactive/length(state.firms.capital)
            row["CES_power_sum"]=direct.capital_power
            row["nearest_threshold_abs_distance"]=
                direct.nearest_threshold_abs_distance
        end
        push!(rows,row)
    end
    if t==5
        count==3 || error("original 20261101 t=5 fixed-point count is no longer 3")
        isempty(positive) ||
            error("unexpected positive-X equilibrium at original benchmark t=5")
        length(check_set)==count ||
            error("not all original t=5 roots were directly checked")
    end
    selected_K=selected==0 ? 0.0 : ws.fixed_capital[selected]
    selected_X=selected==0 ? 0.0 : ws.fixed_investment[selected]
    # Exact identity: Delta K = initial reaction residual + change in Phi.
    reaction_at_expected=fixed.K-expected_K18
    delta_K=selected==0 ? 0.0 : selected_K-expected_K18
    change_in_reaction=selected==0 ? 0.0 : selected_K-fixed.K
    identity_res=selected==0 ? 0.0 :
        delta_K-(reaction_at_expected+change_in_reaction)
    abs(identity_res)<1e-11 || error("exact capital feedback identity failed")
    summary=Dict{String,Any}(
        "t"=>t,"K_prior"=>priorK,"expected_X"=>expected.X_next,
        "K_expected_Eq18"=>expected_K18,"K_expected_AppendixA"=>expected.K_next,
        "n_fixed"=>count,"n_positive_X"=>length(positive),
        "minimum_fixed_X"=>minimum(ws.fixed_investment[1:count]),
        "maximum_fixed_X"=>maximum(ws.fixed_investment[1:count]),
        "positive_X_requires_capital_strictly_above"=>cutoff,
        "maximum_fixed_K_minus_positive_X_cutoff"=>
            ws.fixed_capital[count]-cutoff,
        "nearest_fixed_rank"=>nearest,
        "nearest_positive_fixed_rank"=>closest_positive,
        "selected_rank"=>selected,
        "K_ES1"=>selected==0 ? "undefined" : selected_K,
        "X_ES1"=>selected==0 ? "undefined" : selected_X,
        "X_fixed_at_expected"=>fixed.X,
        "Phi_at_expected"=>fixed.K,
        "reaction_residual_at_expected"=>reaction_at_expected,
        "selected_minus_expected_K"=>selected==0 ? "undefined" : delta_K,
        "Phi_selected_minus_Phi_expected"=>
            selected==0 ? "undefined" : change_in_reaction,
        "capital_feedback_identity_residual"=>identity_res,
        "max_direct_K_residual"=>max_verified_residual,
        "max_direct_X_residual"=>max_verified_x_residual,
        "directly_verified_roots"=>length(check_set),
        "positive_investment_failure"=>isempty(positive),
        "ES1_log_investment_defined"=>!isempty(positive),
        "original_iid_shock_not_rejected"=>true
    )
    return (summary=summary,rows=rows)
end

function stage10_feasibility_at_frozen_state!(state, draws::Int, seed::Int)
    state.period==5 || error("feasibility experiment must be at original t=5")
    state.firms.next_ready || error("original t=5 shock has not been drawn")
    firms=state.firms
    cal=state.calibration
    original_a=copy(firms.next_productivity)
    original_K=state.K
    original_mu=state.mu_prior
    original_capital=copy(firms.capital)
    original_current_a=copy(firms.productivity)
    ws=state.workspace
    orig_result=enumerate_fixed_points!(ws,firms,cal)
    orig_roots=[(ws.fixed_capital[j],ws.fixed_investment[j])
                for j in 1:orig_result.total]
    rng=Xoshiro(UInt64(seed) ⊻ STAGE10_MC_TAG)
    rows=Vector{Dict{String,Any}}()
    feasible=0
    c=cal.constants
    for j in 1:draws
        randn!(rng,firms.next_productivity)
        @inbounds @simd for i in eachindex(firms.next_productivity)
            firms.next_productivity[i]=exp(
                c.mean_log_a+cal.parameters.sigma_log_a*firms.next_productivity[i]
            )
        end
        roots=enumerate_fixed_points!(ws,firms,cal)
        feasible+=(roots.admissible>0)
        push!(rows,Dict{String,Any}(
            "draw"=>j,"feasible"=>roots.admissible>0,
            "n_fixed"=>roots.total,"n_positive_X"=>roots.admissible,
            "max_fixed_X"=>maximum(ws.fixed_investment[1:roots.total]),
            "min_fixed_X"=>minimum(ws.fixed_investment[1:roots.total])
        ))
    end
    copyto!(firms.next_productivity,original_a)
    check=enumerate_fixed_points!(ws,firms,cal)
    rerun=[(ws.fixed_capital[j],ws.fixed_investment[j]) for j in 1:check.total]
    isequal(rerun,orig_roots) ||
        error("original t=5 shock fixed-point correspondence was not restored")
    state.K==original_K && state.mu_prior==original_mu &&
    firms.capital==original_capital &&
    firms.productivity==original_current_a &&
    firms.next_ready ||
        error("conditional iid MC modified a predetermined original benchmark state")
    ci=stage10_wilson(feasible,draws)
    return (summary=Dict{String,Any}(
        "frozen_date"=>5,"draws"=>draws,"feasible_draws"=>feasible,
        "infeasible_draws"=>draws-feasible,
        "p_feasible"=>feasible/draws,
        "MC_SE"=>sqrt((feasible/draws)*(1-feasible/draws)/draws),
        "wilson95_low"=>ci.low,"wilson95_high"=>ci.high,
        "original_t5_shock_untouched_and_restored"=>true,
        "conditioning"=>"Original no-rejection t=5 predetermined firm state; iid draws only for a separate diagnostic, NOT a replacement of original t=5 shock."
    ),rows=rows)
end

function stage10_csv(path::String,rows::Vector{Dict{String,Any}},cols)
    open(path,"w") do io
        println(io,join(cols,","))
        for row in rows
            println(io,join((string(get(row,col,"")) for col in cols),","))
        end
    end
end

function run_stage10()
    n=length(ARGS)>=1 ? parse(Int,ARGS[1]) : 110000
    seed=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 20261101
    draws=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 32
    n>0 && draws>=2 || error("N must be positive and draws >= 2")
    geometry=Vector{Dict{String,Any}}()
    fixed_rows=Vector{Dict{String,Any}}()
    for date in (0,4,5)
        state=stage10_benchmark_state(date,n,seed)
        snap=stage10_geometry!(
            state,date;independently_check_every_root=(date==5)
        )
        append!(fixed_rows,snap.rows)
        push!(geometry,snap.summary)
        println("STAGE10 GEOMETRY date=",date," summary=",snap.summary)
        if date==5
            mc=stage10_feasibility_at_frozen_state!(state,draws,seed)
            println("STAGE10 t5 FEASIBILITY summary=",mc.summary)
            root=normpath(joinpath(@__DIR__,"..","results"))
            mkpath(root)
            stem="diagnostic_equilibrium_feasibility_N$(n)_seed$(seed)_draws$(draws)"
            stage10_csv(joinpath(root,stem*"_fixed_points.csv"),fixed_rows,(
                "t","rank","K","X","positive_X","X_minus_forecast",
                "logK_minus_forecast","nearest_to_expected",
                "nearest_positive_to_expected","ES1_selected",
                "direct_K_residual","direct_X_residual",
                "investment_share","divestment_share","inactive_share",
                "CES_power_sum","nearest_threshold_abs_distance"
            ))
            stage10_csv(joinpath(root,stem*"_feasibility_draws.csv"),mc.rows,(
                "draw","feasible","n_fixed","n_positive_X",
                "max_fixed_X","min_fixed_X"
            ))
            open(joinpath(root,stem*".toml"),"w") do io
                TOML.print(io,Dict(
                    "experiment"=>"issue5_stage10_original_no_rejection_t5",
                    "not_published_replication"=>true,
                    "N"=>n,"seed"=>seed,"draws"=>draws,
                    "julia_version"=>string(VERSION),
                    "geometry"=>geometry,"t5_feasibility_MC"=>mc.summary
                ))
            end
            println("STAGE10 COMPLETED original t5 certificate and frozen iid MC")
            println("STAGE10 OUTPUT ",stem)
        end
    end
    return nothing
end

if abspath(PROGRAM_FILE)==@__FILE__
    run_stage10()
end
