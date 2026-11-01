
"""
Independent post-solve checks are carried out on the final accepted
aggregate path, not on the root solver's termination flag alone.
"""
function transition_diagnostics(sol::TransitionSolveResult{T};
                                compute_euler::Bool=true) where {T}
    sol.converged || throw(ArgumentError("nonconverged solution cannot be audited as equilibrium"))
    w=sol.workspace
    input=sol.input
    n=length(input.B)
    Kerror=zero(T)
    Lerror=zero(T)
    Gerror=zero(T)
    caperror=zero(T)
    maxmass=maximum(w.distribution_error)
    maxbudget=zero(T)
    @inbounds for t in 1:n
        if t<n
            Kerror=max(Kerror,abs((w.aggregate_savings[t]-input.B[t])/w.q[t]-w.K[t+1]))
            investment=w.K[t+1]-(1-input.params.delta)*w.K[t]
            caperror=max(caperror,abs(investment-w.I[t]))
            supply=w.Y[t]-w.I[t]-w.adjustment[t]-
                   input.reference.export_shift*w.terms[t]^(-input.params.theta_star)
        else
            gap=w.I[t-1]/w.K[t]-input.params.delta
            cost=input.phi_ac/2*w.K[t]*gap^2
            supply=w.Y[t]-w.I[t-1]-cost-
                   input.reference.export_shift*w.terms[t]^(-input.params.theta_star)
        end
        Gerror=max(Gerror,abs(supply-w.home_budget[t]))
        Lerror=max(Lerror,abs(w.labor_implied[t]-w.L[t]))
        if t<n
            m=consumption_cost(w.terms[t],input.params)
            gross=1+w.rate[t]
            for j in axes(w.fine_savings,2), i in axes(w.fine_savings,1)
                a=input.grids.fine[i]
                earnings=w.wage[t]*w.L[t]*input.chain.levels[j]+w.profits[t]
                budget=gross*a+earnings-
                       m*w.fine_consumption[i,j,t]-w.fine_savings[i,j,t]
                maxbudget=max(maxbudget,abs(budget))
            end
        end
    end
    euler_max=zero(T)
    euler_weighted=zero(T)
    euler_mass=zero(T)
    if compute_euler
        nf=length(input.grids.fine)
        ns=length(input.chain.levels)
        slopes=Matrix{T}(undef,nf,ns)
        @inbounds for t in 1:n-1
            future=view(w.fine_consumption,:,:,t+1)
            for j in 1:ns
                pchip_slopes!(view(slopes,:,j),input.grids.fine,
                              view(future,:,j))
            end
            current=view(w.fine_consumption,:,:,t)
            decisions=view(w.fine_savings,:,:,t)
            mu=view(w.distribution,:,:,t)
            threshold=input.grids.borrowing_limit+
                      T(1e-8)*max(one(T),abs(input.grids.borrowing_limit))
            for j in 1:ns, i in 1:nf
                c=current[i,j]
                saving=decisions[i,j]
                lhs=c^(-input.params.gamma)
                expected=zero(T)
                for jp in 1:ns
                    cp=pchip_evaluate(input.grids.fine,view(future,:,jp),
                                      view(slopes,:,jp),saving)
                    cp>0 || throw(DomainError(cp,"nonpositive future consumption in Euler audit"))
                    expected+=input.chain.transition[j,jp]*cp^(-input.params.gamma)
                end
                rhs=input.initial.beta*(1+w.rate[t+1])*expected
                gap=(lhs-rhs)/lhs
                violation=saving<=threshold ? max(zero(T),-gap) : abs(gap)
                euler_max=max(euler_max,violation)
                euler_weighted+=mu[i,j]*violation
                euler_mass+=mu[i,j]
            end
        end
    end
    return (max_asset_market_level=Kerror,
            max_labor_market_level=Lerror,
            max_goods_market_level=Gerror,
            max_capital_accumulation=caperror,
            max_household_budget=maxbudget,
            max_probability_mass_deviation=maxmass,
            max_dynamic_euler=euler_max,
            weighted_mean_dynamic_euler=compute_euler ?
                euler_weighted/euler_mass : T(NaN),
            terminal_capital_gap=abs(w.K[end]-input.terminal.K),
            terminal_labor_gap=abs(w.L[end]-input.terminal.L),
            max_scaled_solver_residual=sol.maximum_scaled_residual,
            iterations=sol.iterations,
            evaluations=sol.evaluations)
end

function transition_stationary_check(input::TransitionInputs{T}) where {T}
    n=length(input.B)
    constant=TransitionInputs(input.params,input.reference,
                              input.chain,input.grids,input.initial,input.initial,
                              fill(input.initial.B,n),ones(T,n),
                              input.phi_ac,one(T),input.initial.B)
    x=vcat(fill(input.initial.K,n-1),fill(input.initial.L,n))
    res=zeros(T,2n-1)
    workspace=TransitionWorkspace(constant)
    transition_residual!(res,x,constant,workspace)
    return (max_scaled_residual=maximum(abs,res),
            max_scaled_asset=maximum(abs,view(res,1:n-1)),
            max_scaled_labor=maximum(abs,view(res,n:2n-1)),
            initial_labor=constant.initial.L,
            implied_labor=workspace.labor_implied[1],
            note="Author's transition labor formula omits 1/varphi power; "*
                 "initial steady-state stationary limit can fail even with correct code")
end

function export_transition_csv(sol::TransitionSolveResult{T},
                               filename::AbstractString) where {T}
    sol.converged || throw(ArgumentError("do not export a nonconverged equilibrium"))
    w=sol.workspace
    input=sol.input
    mkpath(dirname(filename))
    open(filename,"w") do io
        println(io,"quarter,B,K,L,Y,I,q,rK,r,w,terms,nominal_exchange_rate,profits,C_H,C_F,C_aggregate,A_supply,goods_residual,capital_residual,labor_residual,distribution_mass")
        n=length(input.B)
        p=input.params
        @inbounds for t in 1:n
            home=w.aggregate_home[t]
            foreign=(1-p.chi)/p.chi*w.terms[t]^p.theta*home
            # Cobb-Douglas aggregate (multiplicative constant immaterial for IRFs).
            consum=home^p.chi*foreign^(1-p.chi)
            goods=if t<n
                w.Y[t]-w.I[t]-w.adjustment[t]-
                input.reference.export_shift*w.terms[t]^(-p.theta_star)-home
            else
                gap=w.I[t-1]/w.K[t]-p.delta
                terminal_cost=input.phi_ac/2*w.K[t]*gap^2
                w.Y[t]-w.I[t-1]-terminal_cost-
                input.reference.export_shift*w.terms[t]^(-p.theta_star)-home
            end
            capital=t<n ? (w.aggregate_savings[t]-input.B[t])/w.q[t]-w.K[t+1] : T(NaN)
            labor=w.labor_implied[t]-w.L[t]
            println(io,join((t-1,input.B[t],w.K[t],w.L[t],w.Y[t],w.I[t],
                             w.q[t],w.rental[t],w.rate[t],w.wage[t],
                             w.terms[t],inv(w.terms[t]),w.profits[t],
                             home,foreign,consum,w.aggregate_savings[t],
                             goods,capital,labor,sum(view(w.distribution,:,:,t))),","))
        end
    end
    return filename
end

# SVG is emitted directly from Julia so plotting has no dependency on Python,
# R, or external libraries, and no unverified Julia plotting API is needed.
function svg_line_panel(io,title::String,ylab::String,values::AbstractVector{T},
                        x0::Float64,y0::Float64,width::Float64,height::Float64;
                        reference::Float64=0.0,firstquarter::Int=0) where {T}
    left=x0+55
    top=y0+36
    w=width-73
    h=height-79
    n=length(values)
    n>=2 || throw(ArgumentError("need at least two chart points"))
    vals=Float64.(values)
    all(isfinite,vals) || throw(ArgumentError("nonfinite plot value"))
    ymin=min(minimum(vals),reference)
    ymax=max(maximum(vals),reference)
    delta=ymax-ymin
    pad=max(0.07*delta,0.25)
    ymin-=pad
    ymax+=pad
    conv_x(i)=left+(i-1)/(n-1)*w
    conv_y(y)=top+(ymax-y)/(ymax-ymin)*h
    println(io,"<rect x='$(x0)' y='$(y0)' width='$(width)' height='$(height)' fill='white'/>")
    println(io,"<text x='$(x0+width/2)' y='$(y0+20)' class='heading' text-anchor='middle'>$(title)</text>")
    println(io,"<path d='M $left $top v $h h $w' fill='none' stroke='#333' stroke-width='1'/>")
    ypos=conv_y(reference)
    println(io,"<line x1='$left' y1='$ypos' x2='$(left+w)' y2='$ypos' stroke='#a3a3a3' stroke-dasharray='4,4'/>")
    print(io,"<polyline fill='none' stroke='#2563eb' stroke-width='2.7' points='")
    for i in 1:n
        print(io,conv_x(i),",",conv_y(vals[i])," ")
    end
    println(io,"'/>")
    for k in 0:4
        xx=left+w*k/4
        tick=firstquarter+round(Int,(n-1)*k/4)
        println(io,"<text x='$xx' y='$(top+h+16)' class='tick' text-anchor='middle'>$tick</text>")
    end
    for k in 0:4
        yy=top+h*k/4
        tick=ymax-(ymax-ymin)*k/4
        number=round(tick;digits=2)
        println(io,"<text x='$(left-7)' y='$(yy+4)' class='tick' text-anchor='end'>$number</text>")
    end
    println(io,"<text x='$(x0+width/2)' y='$(y0+height-6)' class='tick' text-anchor='middle'>Quarters</text>")
    println(io,"<text x='$(left)' y='$(y0+32)' class='tick'>$ylab</text>")
end

function write_transition_svg(sol::TransitionSolveResult{T},destination::AbstractString;
                              horizon::Int=100) where {T}
    sol.converged || throw(ArgumentError("cannot plot an unconverged transition"))
    w=sol.workspace
    input=sol.input
    p=input.params
    mkpath(destination)
    n=min(horizon,length(input.B))
    n>=3 || throw(ArgumentError("horizon too short for SVG"))
    states=collect(1:n-1)
    # Author Figures 2/3 place initial SS at plotted quarter zero,
    # followed by model transition period 1 at plotted quarter one.
    pct(series,initial)=vcat(0.0,
        [100*log(series[t]/initial) for t in states])
    growth(series,initial)=vcat(0.0,
        [100*(series[t]/initial-1) for t in states])
    home=collect(w.aggregate_home)
    foreign=[(1-p.chi)/p.chi*w.terms[t]*home[t] for t in eachindex(home)]
    foreign0=(1-p.chi)/p.chi*input.initial.terms*input.initial.consumption_home
    makefig(file,title,panels)=open(joinpath(destination,file),"w") do io
        println(io,"<svg xmlns='http://www.w3.org/2000/svg' width='1120' height='780' viewBox='0 0 1120 780'>")
        println(io,"<style>.heading{font: 600 18px 'Cambria Math',Cambria,serif;} .tick{font: 13px 'Cambria Math',Cambria,serif;} .caption{font: 17px Cambria,serif;}</style>")
        println(io,"<rect width='1120' height='780' fill='white'/>")
        println(io,"<text x='560' y='29' text-anchor='middle' class='heading'>$title</text>")
        for (i,(name,units,data,baseline)) in enumerate(panels)
            row=(i-1)÷2
            col=(i-1)%2
            svg_line_panel(io,name,units,data,20.0+col*550.0,
                           49.0+row*355.0,535.0,343.0;reference=baseline)
        end
        println(io,"<text x='560' y='769' class='caption' text-anchor='middle'>de Ferra, Mitman &amp; Romei (2020) — Julia replication; source: solved nonlinear transition</text>")
        println(io,"</svg>")
    end
    # Figure 1 in the paper contains ALSO an unanticipated contraction
    # 40 quarters later. Stage 2 plots expansion only; label it as partial.
    open(joinpath(destination,"figure_1.svg"),"w") do io
        println(io,"<svg xmlns='http://www.w3.org/2000/svg' width='850' height='500' viewBox='0 0 850 500'>")
        println(io,"<style>.heading{font: 600 23px Cambria,serif;} .tick{font: 16px Cambria,serif;}</style>")
        println(io,"<rect width='850' height='500' fill='white'/>")
        println(io,"<text x='425' y='30' class='heading' text-anchor='middle'>Foreign Credit Supply — expansion only</text>")
        expansion=[-0.25*input.B[t]/input.initial.Y for t in 1:n]
        svg_line_panel(io,"Figure 1: Stage 2 component","Fraction of annual GDP",
                       expansion,35.0,56.0,775.0,380.0;reference=0.0)
        println(io,"<text x='425' y='477' class='tick' text-anchor='middle'>Unexpected contraction is excluded from Stage 2 and will be added in Stage 3.</text>")
        println(io,"</svg>")
    end
    makefig("figure_2.svg","Consumption, output, and investment — credit expansion",[
        ("Home-good consumption Cₕ","100 × log deviation (%)",
          pct(home,input.initial.consumption_home),0.0),
        ("Foreign-good consumption C𝒻","100 × log deviation (%)",
          pct(foreign,foreign0),0.0),
        ("Output Y","100 × log deviation (%)",pct(w.Y,input.initial.Y),0.0),
        ("Investment I","Percent from initial SS",
          growth(w.I,input.initial.investment),0.0)
    ])
    makefig("figure_3.svg","Trade prices, real wage, exchange rate, and interest",[
        ("Terms of trade Pₕ/P𝒻","100 × log deviation (%)",
          pct(w.terms,input.initial.terms),0.0),
        ("Real wage w","100 × log deviation (%)",
          pct(w.wage,input.initial.wage),0.0),
        ("Nominal exchange rate 1/x","100 × log deviation (%)",
          pct(inv.(w.terms),inv(input.initial.terms)),0.0),
        ("Real interest rate r","Quarterly percentage rate",
          vcat(100*input.initial.r,[100*w.rate[t] for t in states]),
          100*input.initial.r)
    ])
    return (joinpath(destination,"figure_1.svg"),
            joinpath(destination,"figure_2.svg"),
            joinpath(destination,"figure_3.svg"))
end
