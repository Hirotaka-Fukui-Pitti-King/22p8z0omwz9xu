#!/usr/bin/env julia
# Nirei (2015), Issue #5, Stage 12: DIAGNOSTIC choice-set extension only.
# This file does NOT modify the original three-choice production equilibrium.
#
# julia --project=. scripts/diagnostic_multistep_equilibrium.jl small
# julia --project=. scripts/diagnostic_multistep_equilibrium.jl benchmark [N] [seed]
#
# Exact interval enumeration of K = Phi_J(K), for kappa in -J:J.
# Distinguish an admissible finite-J choice-set counterfactual from the
# published one-notch restriction and from a full nonlinear RE equilibrium.
include(joinpath(@__DIR__, "diagnostic_equilibrium_feasibility.jl"))

struct MultistepRoot
    K::Float64
    X::Float64
end

struct MultistepEnumeration
    J::Int
    N::Int
    roots::Vector{MultistepRoot}
    switches::Vector{Float64}
    order::Vector{Int}
end

@inline function multistep_policy_at_K(
    K::Float64, d::Float64, lambda_i::Float64, b::Float64, J::Int
)
    # Evaluate the breakpoint inequalities rather than taking ceil(log(.)).
    # At a threshold the OLD policy applies, exactly as in Eq. (16).
    exponent=-J
    @inbounds for j in -J:(J-1)
        K > d*lambda_i^j/b && (exponent+=1)
    end
    return exponent
end

function multistep_direct_response(firms, calibration, K::Float64, J::Int)
    J>=1 || error("J must be >= 1")
    K>0 && isfinite(K) || error("candidate capital must be finite and positive")
    n=length(firms.capital)
    p=calibration.parameters
    rho=calibration.constants.rho
    ces=0.0
    total_x=0.0
    investing_one=0
    investing_many=0
    divesting_one=0
    divesting_many=0
    inactive=0
    nearest=Inf
    @inbounds for i in 1:n
        d=(1.0-p.delta)*firms.capital[i]
        a=firms.next_productivity[i]
        li=firms.lambda[i]
        b=threshold_coefficient(li,a,calibration)
        exponent=multistep_policy_at_K(K,d,li,b,J)
        newk=d*li^exponent
        ces+=a^(rho/p.alpha)*newk^rho/n
        total_x+=(newk-d)/n
        investing_one+=(exponent==1)
        investing_many+=(exponent>=2)
        divesting_one+=(exponent==-1)
        divesting_many+=(exponent<=-2)
        inactive+=(exponent==0)
        for j in -J:(J-1)
            nearest=min(nearest,abs(K-d*li^j/b))
        end
    end
    investing_one+investing_many+divesting_one+divesting_many+inactive==n ||
        error("multi-step policy counts do not sum to N")
    return (
        K=ces^(inv(rho)), X=total_x, ces_power=ces,
        invest_one=investing_one, invest_many=investing_many,
        divest_one=divesting_one, divest_many=divesting_many,
        inactive=inactive, nearest_switch=nearest
    )
end

function multistep_enumerate(firms, calibration, J::Int)
    J>=1 || error("J must be >= 1")
    firms.next_ready || error("future productivity has not been drawn")
    n=length(firms.capital)
    m=2J*n
    switches=Vector{Float64}(undef,m)
    dvec=Vector{Float64}(undef,n)
    basevec=Vector{Float64}(undef,n)
    xsum=0.0
    psum=0.0
    rho=calibration.constants.rho
    p=calibration.parameters
    invn=inv(Float64(n))
    @inbounds for i in 1:n
        d=(1-p.delta)*firms.capital[i]
        a=firms.next_productivity[i]
        li=firms.lambda[i]
        b=threshold_coefficient(li,a,calibration)
        (d>0 && a>0 && li>1 && isfinite(b) && b>0) ||
            error("invalid predetermined firm state at index $i")
        base=a^(rho/p.alpha)*d^rho*invn
        dvec[i]=d
        basevec[i]=base
        xsum+=d*(li^(-J)-1)*invn
        psum+=base*li^(-J*rho)
        for notch in -J:(J-1)
            pos=(i-1)*(2J)+notch+J+1
            val=d*li^notch/b
            val>0 && isfinite(val) || error("invalid breakpoint at firm $i")
            switches[pos]=val
        end
    end
    order=sortperm(switches)
    roots=MultistepRoot[]
    left=0.0
    pos=1
    @inbounds while pos<=m
        boundary=switches[order[pos]]
        candidate=psum^(inv(rho))
        if left<candidate<=boundary
            push!(roots,MultistepRoot(candidate,xsum))
        end
        while pos<=m && switches[order[pos]]==boundary
            event=order[pos]
            i=div(event-1,2J)+1
            step=mod(event-1,2J)+1
            notch=step-J-1
            li=firms.lambda[i]
            xsum+=dvec[i]*li^notch*(li-1)*invn
            psum+=basevec[i]*li^(rho*notch)*(li^rho-1)
            pos+=1
        end
        left=boundary
    end
    candidate=psum^(inv(rho))
    if left<candidate
        push!(roots,MultistepRoot(candidate,xsum))
    end
    isempty(roots) && error("no multi-step fixed point found (J=$J)")
    return MultistepEnumeration(J,n,roots,switches,order)
end

function multistep_validate_roots!(enumeration, firms, cal;
                                   all_roots::Bool=true)
    roots=enumeration.roots
    indices=all_roots ? collect(eachindex(roots)) :
        unique([firstindex(roots),lastindex(roots)])
    max_k=0.0
    max_x=0.0
    for idx in indices
        root=roots[idx]
        direct=multistep_direct_response(firms,cal,root.K,enumeration.J)
        errk=abs(direct.K-root.K)
        errx=abs(direct.X-root.X)
        max_k=max(max_k,errk)
        max_x=max(max_x,errx)
        isapprox(root.K,direct.K;atol=2e-9,rtol=1e-10) ||
            error("direct capital residual at J=$(enumeration.J), root $idx: $errk")
        isapprox(root.X,direct.X;atol=2e-9,rtol=1e-10) ||
            error("direct investment residual at J=$(enumeration.J), root $idx: $errx")
    end
    return (n_checked=length(indices),max_k_residual=max_k,max_x_residual=max_x)
end

# Independent O(N*(2JN)) oracle is used ONLY with N <= 50.
# The oracle forms candidate Phi from direct firm choices on each
# breakpoint interval; it does not reuse the cumulative jump accounting.
function multistep_small_oracle(firms,cal,J::Int)
    n=length(firms.capital)
    n<=50 || error("independent brute force oracle must use N <= 50")
    limits=Float64[]
    for i in 1:n
        d=(1-cal.parameters.delta)*firms.capital[i]
        li=firms.lambda[i]
        b=threshold_coefficient(li,firms.next_productivity[i],cal)
        for j in -J:(J-1)
            push!(limits,d*li^j/b)
        end
    end
    sort!(limits)
    unique!(limits)
    oracle=MultistepRoot[]
    left=0.0
    for right in vcat(limits,Inf)
        probe=isfinite(right) ?
            (left==0 ? right/2 : left+(right-left)/2) :
            left+max(1.0,left/2)
        response=multistep_direct_response(firms,cal,probe,J)
        candidate=response.K
        if left<candidate<=right
            atroot=multistep_direct_response(firms,cal,candidate,J)
            isapprox(atroot.K,candidate;atol=2e-9,rtol=1e-10) ||
                error("small oracle candidate is not a fixed point")
            push!(oracle,MultistepRoot(candidate,atroot.X))
        end
        left=right
    end
    return oracle
end

function stage12_small_tests()
    for seed in (20261101,20261102,20261103,20261104)
        st=initialize_simulation(;N=24,base_seed=seed)
        draw_productivity!(st.firms,st.rngs,st.calibration)
        old=enumerate_fixed_points!(st.workspace,st.firms,st.calibration)
        oldks=copy(st.workspace.fixed_capital[1:old.total])
        oldxs=copy(st.workspace.fixed_investment[1:old.total])
        for J in 1:4
            new=multistep_enumerate(st.firms,st.calibration,J)
            independent=multistep_small_oracle(st.firms,st.calibration,J)
            length(new.roots)==length(independent) ||
                error("independent oracle root count differs for seed=$seed J=$J")
            for (a,b) in zip(new.roots,independent)
                isapprox(a.K,b.K;atol=2e-9,rtol=1e-10) ||
                    error("small oracle root capital differs for seed=$seed J=$J")
                isapprox(a.X,b.X;atol=2e-9,rtol=1e-10) ||
                    error("small oracle investment differs for seed=$seed J=$J")
            end
            verified=multistep_validate_roots!(new,st.firms,st.calibration)
            if J==1
                length(new.roots)==old.total ||
                    error("J=1 fails production root count regression")
                for (i,root) in enumerate(new.roots)
                    isapprox(root.K,oldks[i];atol=2e-9,rtol=1e-10) ||
                        error("J=1 production capital regression failed")
                    isapprox(root.X,oldxs[i];atol=2e-9,rtol=1e-10) ||
                        error("J=1 production investment regression failed")
                end
            end
            println("STAGE12 SMALL seed=",seed," J=",J,
                    " root_count=",length(new.roots),
                    " direct_max_K=",verified.max_k_residual,
                    " direct_max_X=",verified.max_x_residual)
        end
    end
    println("STAGE12 SMALL ORACLE SUCCESS")
end

function stage12_csv(path,rows,cols)
    open(path,"w") do io
        println(io,join(cols,","))
        for row in rows
            println(io,join((string(get(row,col,"")) for col in cols),","))
        end
    end
end

function stage12_benchmark(N::Int,seed::Int)
    N>0 || error("N must be positive")
    st=stage10_benchmark_state(5,N,seed)
    draw_productivity!(st.firms,st.rngs,st.calibration)
    expected=compute_expected_equilibrium(
        st.re,st.calibration,st.K,st.mu_prior
    )
    ce=st.calibration
    ke=(1-ce.parameters.delta)*st.K+
       ce.threshold.A^(ce.constants.rho-1)/ce.constants.rho*expected.X_next
    # Save original production one-notch roots BEFORE the extended sweeps.
    legacy=enumerate_fixed_points!(st.workspace,st.firms,ce)
    legacyK=copy(st.workspace.fixed_capital[1:legacy.total])
    legacyX=copy(st.workspace.fixed_investment[1:legacy.total])
    if N==110000 && seed==20261101
        legacy.total==3 && legacy.admissible==0 ||
            error("Original benchmark t5 root count/feasibility regression failed")
        isapprox(maximum(legacyX),-0.057398135458978324;atol=5e-10) ||
            error("original benchmark t5 max investment regression failed")
    end
    summaries=Dict{String,Any}[]
    rows=Dict{String,Any}[]
    for J in 1:4
        out=multistep_enumerate(st.firms,ce,J)
        verification=multistep_validate_roots!(out,st.firms,ce;all_roots=true)
        if J==1
            length(out.roots)==legacy.total ||
                error("benchmark J1 original fixed point count mismatch")
            for (i,root) in enumerate(out.roots)
                isapprox(root.K,legacyK[i];atol=2e-9,rtol=1e-10) ||
                    error("benchmark J1 original capital mismatch")
                isapprox(root.X,legacyX[i];atol=2e-9,rtol=1e-10) ||
                    error("benchmark J1 original investment mismatch")
            end
        end
        feasible=[i for i in eachindex(out.roots) if out.roots[i].X>0.0]
        selected=0
        score=Inf
        for i in feasible
            crit=abs(log(out.roots[i].X/expected.X_next))
            if crit<score
                selected=i
                score=crit
            end
        end
        maxX=maximum(root.X for root in out.roots)
        near=argmin([abs(log(root.K/ke)) for root in out.roots])
        representative=selected>0 ? selected : length(out.roots)
        selected_response=multistep_direct_response(
            st.firms,ce,out.roots[representative].K,J
        )
        info=Dict{String,Any}(
            "J"=>J,"N"=>N,"t"=>5,"seed"=>seed,
            "n_fixed"=>length(out.roots),"n_positive_X"=>length(feasible),
            "maximum_equilibrium_X"=>maxX,
            "minimum_equilibrium_X"=>minimum(root.X for root in out.roots),
            "minimum_equilibrium_K"=>first(out.roots).K,
            "maximum_equilibrium_K"=>last(out.roots).K,
            "nearest_expected_K_rank"=>near,
            "nearest_expected_K"=>out.roots[near].K,
            "comparison_ES1_rank"=>selected,
            "comparison_ES1_K"=>selected==0 ? "undefined" : out.roots[selected].K,
            "comparison_ES1_X"=>selected==0 ? "undefined" : out.roots[selected].X,
            "representative_root_rank"=>representative,
            "representative_invest_one"=>selected_response.invest_one,
            "representative_invest_2plus"=>selected_response.invest_many,
            "representative_divest_one"=>selected_response.divest_one,
            "representative_divest_2plus"=>selected_response.divest_many,
            "representative_no_adjustment"=>selected_response.inactive,
            "expected_X"=>expected.X_next,"expected_K_Eq18"=>ke,
            "all_roots_independently_checked"=>verification.n_checked,
            "max_direct_K_residual"=>verification.max_k_residual,
            "max_direct_X_residual"=>verification.max_x_residual,
            "diagnostic_only"=>true
        )
        for (idx,root) in enumerate(out.roots)
            push!(rows,Dict{String,Any}(
                "J"=>J,"rank"=>idx,"K"=>root.K,"X"=>root.X,
                "positive_X"=>root.X>0,"comparison_ES1"=>idx==selected,
                "near_expected_K"=>idx==near
            ))
        end
        push!(summaries,info)
        println("STAGE12 BENCHMARK ",info)
    end
    root=normpath(joinpath(@__DIR__,"..","results"))
    mkpath(root)
    stem="diagnostic_multistep_equilibrium_N$(N)_seed$(seed)"
    stage12_csv(joinpath(root,stem*"_fixed_points.csv"),rows,(
        "J","rank","K","X","positive_X","comparison_ES1","near_expected_K"
    ))
    open(joinpath(root,stem*".toml"),"w") do io
        TOML.print(io,Dict(
            "experiment"=>"issue5_stage12_multistep_choice_set_counterfactual",
            "original_t5_shock_preserved"=>true,
            "not_original_three_choice_replication_for_J_above_1"=>true,
            "N"=>N,"seed"=>seed,
            "julia_version"=>string(VERSION),"summaries"=>summaries
        ))
    end
    println("STAGE12 BENCHMARK COMPLETED output=",stem)
end

if abspath(PROGRAM_FILE)==@__FILE__
    mode=isempty(ARGS) ? "small" : ARGS[1]
    if mode=="small"
        stage12_small_tests()
    elseif mode=="benchmark"
        n=length(ARGS)>=2 ? parse(Int,ARGS[2]) : 110000
        seed=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 20261101
        stage12_benchmark(n,seed)
    else
        error("Usage: diagnostic_multistep_equilibrium.jl small | benchmark [N] [seed]")
    end
end
