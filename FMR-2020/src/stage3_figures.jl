
"""
Original Figures.m panels 1, 4, 5: two-source paths (pre-shock boom and
new surprise reversal) and the original 30-quarter pre-shock normalization.
All SVG art is generated in Julia, without Python/MATLAB/R dependencies.
"""
function stage3_svg_open(io,caption::String,width::Int,height::Int)
    println(io,"<svg xmlns='http://www.w3.org/2000/svg' width='$width' height='$height' viewBox='0 0 $width $height'>")
    println(io,"<style>.title{font:600 19px Cambria,serif}.axis{font:13px Cambria,serif}.legend{font:15px Cambria,serif}</style>")
    println(io,"<rect width='$width' height='$height' fill='white'/>")
    println(io,"<text class='title' x='$(width/2)' y='27' text-anchor='middle'>$caption</text>")
end

function stage3_svg_panel(io,name::String,units::String,
                          series_a::AbstractVector,series_b::AbstractVector,
                          left::Real,top::Real,width::Real,height::Real;
                          show_legend::Bool=false)
    length(series_a)==length(series_b) ||
        throw(DimensionMismatch("comparison figure paths"))
    n=length(series_a)
    n>=2 || throw(ArgumentError("plot horizon must be at least two"))
    all(isfinite,series_a) && all(isfinite,series_b) ||
        throw(ArgumentError("SVG cannot plot nonfinite nonlinear GE paths"))
    lo=min(minimum(series_a),minimum(series_b),0.0)
    hi=max(maximum(series_a),maximum(series_b),0.0)
    d=hi-lo
    pad=max(0.12*d,0.4)
    lo-=pad; hi+=pad
    x0=left+59
    y0=top+34
    ww=width-80
    hh=height-83
    xp(i)=x0+(i-1)*ww/(n-1)
    yp(v)=y0+(hi-v)*hh/(hi-lo)
    println(io,"<text class='title' x='$(left+width/2)' y='$(top+20)' text-anchor='middle'>$name</text>")
    println(io,"<path d='M $x0 $y0 v $hh h $ww' fill='none' stroke='#555'/>")
    println(io,"<line x1='$x0' x2='$(x0+ww)' y1='$(yp(0))' y2='$(yp(0))' stroke='#999' stroke-dasharray='3 5'/>")
    for (series,color,dash) in ((series_a,"#2563eb",""),(series_b,"#dc2626","6 4"))
        pts=join((string(xp(i),",",yp(series[i])) for i in 1:n)," ")
        style=dash=="" ? "" : "stroke-dasharray='$dash'"
        println(io,"<polyline fill='none' stroke='$color' stroke-width='2.5' $style points='$pts'/>")
    end
    for k in 0:3
        xx=x0+k*ww/3
        tt=round(Int,k*(n-1)/3)
        yy=y0+k*hh/3
        val=round(hi-(hi-lo)*k/3;digits=2)
        println(io,"<text class='axis' x='$xx' y='$(y0+hh+20)' text-anchor='middle'>$tt</text>")
        println(io,"<text class='axis' x='$(x0-7)' y='$(yy+4)' text-anchor='end'>$val</text>")
    end
    println(io,"<text class='axis' x='$x0' y='$(top+32)'>$units</text>")
    if show_legend
        y=top+height-22
        println(io,"<line x1='$(x0+90)' x2='$(x0+124)' y1='$y' y2='$y' stroke='#2563eb' stroke-width='3'/>")
        println(io,"<text class='legend' x='$(x0+130)' y='$(y+5)'>Flexible</text>")
        println(io,"<line x1='$(x0+260)' x2='$(x0+294)' y1='$y' y2='$y' stroke='#dc2626' stroke-width='3' stroke-dasharray='6 4'/>")
        println(io,"<text class='legend' x='$(x0+300)' y='$(y+5)'>Fixed</text>")
    end
end

function stage3_figure_1(boom::TransitionSolveResult{T},
                          reversal::AbstractVector{T},
                          filename::AbstractString;shock_start::Int=41) where {T}
    n=length(boom.input.B)
    length(reversal)>=n || throw(DimensionMismatch("reversal credit supply"))
    y=boom.input.initial.Y
    supply=-T(0.25).*boom.input.B ./y
    surprise=Vector{T}(undef,n)
    @inbounds for t in 1:n
        surprise[t]=t<shock_start ? supply[t] :
                     -T(0.25)*reversal[t-shock_start+1]/y
    end
    mkpath(dirname(filename))
    open(filename,"w") do io
        stage3_svg_open(io,"Figure 1 — Foreign credit: expansion and unexpected reversal",900,510)
        stage3_svg_panel(io,"Foreign credit supply","Fraction of annual GDP",
                         supply,surprise,35,48,820,395)
        xx=35+59+(shock_start-1)/(n-1)*(820-80)
        println(io,"<line x1='$xx' x2='$xx' y1='82' y2='395' stroke='#e11d48' stroke-dasharray='3 5'/>")
        println(io,"<text class='legend' x='450' y='477' text-anchor='middle'>Blue: initial expansion. Red: unexpected contraction at quarter $(shock_start-1).</text>")
        println(io,"</svg>")
    end
    return filename
end

function stage3_figures_4_5(boom::TransitionSolveResult{T},
                              flexible::TransitionSolveResult{T},
                              fixed::TransitionSolveResult{T},
                              dest::AbstractString;
                              shock_start::Int=41,
                              horizon::Int=30) where {T}
    horizon<=length(flexible.input.B) &&
        horizon<=length(fixed.input.B) ||
        throw(ArgumentError("plot horizon exceeds shock transition"))
    horizon>=3 || throw(ArgumentError("too short IRFs"))
    mkpath(dest)
    pre=boom.workspace
    a=flexible.workspace
    b=fixed.workspace
    p=boom.input.params
    # MATLAB Figures.m shock_start_t=41 refers to the original source
    # 41st transition entry and assumes period zero is that pre-shock
    # anticipated, no-reversal equilibrium value.
    cpref=pre.home_budget[shock_start]
    cfpref=(one(T)-p.chi)/p.chi*
             pre.terms[shock_start]*pre.aggregate_home[shock_start]
    ypref=pre.Y[shock_start]
    prevK=pre.K[shock_start-1]
    invpref=pre.K[shock_start]-(one(T)-p.delta)*prevK
    xpref=pre.terms[shock_start]
    wpref=pre.wage[shock_start]
    rpref=pre.rate[shock_start+1]
    first=1:horizon-1
    logdiff(series,base)=vcat(zero(T),
                             [T(100)*log(series[t]/base) for t in first])
    pctdiff(series,base)=vcat(zero(T),
                             [T(100)*(series[t]/base-one(T)) for t in first])
    foreign_c(w)=[(one(T)-p.chi)/p.chi*w.terms[t]*
                   w.aggregate_home[t] for t in 1:horizon-1]
    fig4=[
        ("Home-good consumption","Log percentage",
         logdiff(a.aggregate_home,cpref),logdiff(b.aggregate_home,cpref)),
        ("Foreign-good consumption","Log percentage",
         logdiff(foreign_c(a),cfpref),logdiff(foreign_c(b),cfpref)),
        ("Output","Log percentage",logdiff(a.Y,ypref),logdiff(b.Y,ypref)),
        ("Investment","Percentage",
         pctdiff(a.I,invpref),pctdiff(b.I,invpref))]
    # Figure 5: the nominal exchange rate under fixed exchange rates
    # follows the model's E_t=E_{-1} restriction. The ToT itself
    # adjusts as domestic prices move; this is NOT a constant-ToT rule.
    e_flex=inv.(a.terms)
    e_initial=inv(xpref)
    e_fixed=fill(e_initial,length(b.terms))
    fig5=[
        ("Terms of trade","Log percentage",
         logdiff(a.terms,xpref),logdiff(b.terms,xpref)),
        ("Real wage","Log percentage",
         logdiff(a.wage,wpref),logdiff(b.wage,wpref)),
        ("Nominal exchange rate","Log percentage",
         logdiff(e_flex,e_initial),logdiff(e_fixed,e_initial)),
        ("Real interest rate","Percentage points",
         vcat(zero(T),[T(100)*(a.rate[t+1]-rpref) for t in first]),
         vcat(zero(T),[T(100)*(b.rate[t+1]-rpref) for t in first]))]
    function draw(file,title,panels)
        destfile=joinpath(dest,file)
        open(destfile,"w") do io
            stage3_svg_open(io,title,1120,800)
            for k in 1:4
                row=(k-1)÷2
                col=(k-1)%2
                name,unit,left,right=panels[k]
                stage3_svg_panel(io,name,unit,left,right,
                                  20+col*550,52+row*357,535,350;
                                  show_legend=k==4)
            end
            println(io,"<text class='axis' x='560' y='784' text-anchor='middle'>Source: de Ferra, Mitman &amp; Romei (2020), Figures.m; shock quarter=41</text>")
            println(io,"</svg>")
        end
        return destfile
    end
    return (draw("figure_4.svg","Figure 4 — Unexpected credit reversal: Flexible vs Fixed",fig4),
            draw("figure_5.svg","Figure 5 — Relative prices, exchange and interest",fig5))
end

function stage3_write_csv(result::TransitionSolveResult{T},
                           ctx::Stage3Context{T},
                           filename::AbstractString;
                           fixed::Bool=false) where {T}
    result.converged || throw(ArgumentError("cannot export nonconverged path"))
    mkpath(dirname(filename))
    w=result.workspace; inp=ctx.economy
    open(filename,"w") do io
        println(io,"period,B,K,L,Y,I,q,r,w,terms,foreign_consumption,home_consumption,wealth,nominal_exchange_rate")
        for t in eachindex(inp.B)
            home=w.aggregate_home[t]
            foreign=(one(T)-inp.params.chi)/inp.params.chi*
                     w.terms[t]*home
            exchange=fixed ? one(T) : inv(w.terms[t])
            println(io,join((t,inp.B[t],w.K[t],w.L[t],w.Y[t],w.I[t],
                w.q[t],w.rate[t],w.wage[t],w.terms[t],foreign,home,
                w.aggregate_savings[t],exchange),","))
        end
    end
    return filename
end
