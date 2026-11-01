# Shape-preserving piecewise cubic Hermite interpolant (Fritsch-Butland
# interior slopes with one-sided limited endpoint derivatives).
# Designed for strictly increasing, nonuniform Julia asset grids.
function pchip_slopes!(d::AbstractVector{T}, x::AbstractVector{T},
                       y::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(x)
    n == length(y) == length(d) || throw(DimensionMismatch("PCHIP input sizes"))
    n >= 2 || throw(ArgumentError("need at least two knots"))
    if n == 2
        d[1] = d[2] = (y[2]-y[1])/(x[2]-x[1])
        return d
    end
    @inbounds for i in 2:(n-1)
        h0 = x[i]-x[i-1]
        h1 = x[i+1]-x[i]
        h0>0 && h1>0 || throw(ArgumentError("nonincreasing interpolation knots"))
        s0 = (y[i]-y[i-1])/h0
        s1 = (y[i+1]-y[i])/h1
        if s0*s1 <= 0
            d[i] = zero(T)
        else
            w0 = 2h1+h0
            w1 = h1+2h0
            d[i] = (w0+w1)/(w0/s0+w1/s1)
        end
    end
    @inbounds begin
        h0 = x[2]-x[1]
        h1 = x[3]-x[2]
        s0 = (y[2]-y[1])/h0
        s1 = (y[3]-y[2])/h1
        endpoint = ((2h0+h1)*s0-h0*s1)/(h0+h1)
        d[1] = sign(endpoint)!=sign(s0) ? zero(T) :
            (sign(s0)!=sign(s1) && abs(endpoint)>3abs(s0) ? 3s0 : endpoint)

        h0 = x[n]-x[n-1]
        h1 = x[n-1]-x[n-2]
        s0 = (y[n]-y[n-1])/h0
        s1 = (y[n-1]-y[n-2])/h1
        endpoint = ((2h0+h1)*s0-h0*s1)/(h0+h1)
        d[n] = sign(endpoint)!=sign(s0) ? zero(T) :
            (sign(s0)!=sign(s1) && abs(endpoint)>3abs(s0) ? 3s0 : endpoint)
    end
    return d
end

@inline function pchip_evaluate(x::AbstractVector{T}, y::AbstractVector{T},
                                d::AbstractVector{T}, z::T) where {T<:AbstractFloat}
    n = length(x)
    if z <= x[1]
        return y[1]+d[1]*(z-x[1])
    elseif z >= x[n]
        return y[n]+d[n]*(z-x[n])
    end
    k = searchsortedlast(x,z)
    @inbounds begin
        h = x[k+1]-x[k]
        u = (z-x[k])/h
        h00 = (1+2u)*(1-u)^2
        h10 = u*(1-u)^2
        h01 = u^2*(3-2u)
        h11 = u^2*(u-1)
        return h00*y[k]+h*h10*d[k]+h01*y[k+1]+h*h11*d[k+1]
    end
end
