
"""
Return the 2T-1-dimensional *economic* perfect-foresight residual:
(T-1) asset-market equations, (T-1) labor-union optimality equations,
and the terminal stationary-labor continuation condition L_T=L_ss.

Paper (2020) eq. (5) and the CRRA/Cobb-Douglas calibration imply
L_t^varphi = chi*(epsilon_w-1)/epsilon_w *
             w_t/(chi_dis*C_H,t).
The author's MATLAB transition erroneously omits the 1/varphi power
and copies C_H,T-1 into the final labor residual. Neither shortcut is
used in this economic closure. Convergence to the terminal steady state
must still be checked independently as the horizon is increased.
"""
function transition_residual!(res::AbstractVector{T},x::AbstractVector{T},
                              input::TransitionInputs{T},
                              w::TransitionWorkspace{T}) where {T}
    n=length(input.B)
    length(res)==2n-1 || throw(DimensionMismatch("2T-1 residual required"))
    w.function_evals+=1
    transition_prices!(w,input,x)
    transition_household_backward!(w,input)
    transition_distribution_forward!(w,input)
    p=input.params
    scaleK=max(one(T),input.initial.K)
    scaleL=max(one(T),input.initial.L)
    @inbounds for t in 1:n-1
        residualK=(w.aggregate_savings[t]-input.B[t])/w.q[t]-w.K[t+1]
        res[t]=residualK/scaleK
    end
    union_factor=p.chi*(p.epsilon_w-1)/p.epsilon_w
    @inbounds for t in 1:n-1
        c=w.aggregate_home[t]
        c>0 || throw(DomainError(c,"nonpositive aggregate home consumption"))
        ratio=union_factor*w.wage[t]/(input.reference.chi_dis*c)
        ratio>zero(T) || throw(DomainError(ratio,"nonpositive union labor FOC"))
        w.labor_implied[t]=ratio^(inv(p.varphi))
        res[n-1+t]=(w.labor_implied[t]-w.L[t])/scaleL
    end
    # The last date is the calibrated stationary continuation, not an
    # additional labor FOC using an arbitrary copy of consumption T-1.
    w.labor_implied[n]=w.L[n]
    res[2n-1]=(w.L[n]-input.terminal.L)/scaleL
    all(isfinite,res) || throw(ArgumentError("nonfinite GE residual"))
    return res
end

function evaluate_transition!(input::TransitionInputs{T}, x::AbstractVector{T};
                              workspace::TransitionWorkspace{T}=TransitionWorkspace(input)) where {T}
    r=zeros(T,length(x))
    transition_residual!(r,x,input,workspace)
    return (residual=r,workspace=workspace)
end

struct TransitionSolveResult{T<:AbstractFloat}
    input::TransitionInputs{T}
    workspace::TransitionWorkspace{T}
    x::Vector{T}
    residual::Vector{T}
    iterations::Int
    evaluations::Int
    maximum_scaled_residual::T
    converged::Bool
end

"""
Dense damped quasi-Newton (Broyden rank-one) with forward finite-difference
Jacobian initialization and periodic refresh. This solves the SAME nonlinear
system as the authors' fsolve; no household approximation, representative-
agent reduction, or frozen distribution is used.

Every trial re-solves prices, backward EGM and forward distribution.
An invalid trial is rejected by step reduction; no economic state is
silently clipped. The caller controls maximum iterations.
"""
function solve_flexible_transition(input::TransitionInputs{T};
                                   x0::AbstractVector{T}=initial_transition_guess(input),
                                   tolerance::T=T(1e-6),
                                   maxiter::Int=70,
                                   refresh_every::Int=15,
                                   verbose::Bool=true) where {T}
    n=2*length(input.B)-1
    length(x0)==n || throw(DimensionMismatch("2T-1 starting vector required"))
    x=Vector{T}(x0)
    w=TransitionWorkspace(input)
    f=zeros(T,n)
    transition_residual!(f,x,input,w)
    ftrial=similar(f)
    xt=similar(x)
    J=Matrix{T}(undef,n,n)
    gradient=zeros(T,n)
    iterations=0
    function jacobian!()
        @inbounds for k in 1:n
            copyto!(xt,x)
            step=T(1e-5)*max(one(T),abs(x[k]))
            xt[k]+=step
            transition_residual!(ftrial,xt,input,w)
            for j in 1:n
                J[j,k]=(ftrial[j]-f[j])/step
            end
        end
    end
    for it in 1:maxiter
        error_max=maximum(abs,f)
        verbose && println("TRANSITION iteration=",it-1,
                           " scaled_max=",error_max,
                           " residual_evals=",w.function_evals)
        if error_max<=tolerance
            # Recompute workspace at accepted candidate after Jacobian
            # probes (which have overwritten mutable state).
            transition_residual!(f,x,input,w)
            return TransitionSolveResult(input,w,x,f,it-1,
                                         w.function_evals,maximum(abs,f),true)
        end
        if it==1 || (it-1)%refresh_every==0
            jacobian!()
            # restore evaluation state: J probes mutated the workspace
            transition_residual!(f,x,input,w)
        end
        gradient .= -(J\f)
        all(isfinite,gradient) ||
            throw(ErrorException("nonfinite GE Newton direction"))
        accepted=false
        alpha=one(T)
        previous_norm=norm(f)
        for _ in 1:16
            @inbounds for k in 1:n
                xt[k]=x[k]+alpha*gradient[k]
            end
            good=true
            try
                transition_residual!(ftrial,xt,input,w)
            catch err
                if err isa DomainError || err isa ArgumentError
                    good=false
                else
                    rethrow()
                end
            end
            if good && norm(ftrial)<(one(T)-T(1e-4)*alpha)*previous_norm
                accepted=true
                break
            end
            alpha/=2
        end
        accepted || throw(ErrorException(
            "nonlinear line search failed at iteration $it; norm=$(previous_norm)"))
        # Rank-one direct-Jacobian Broyden update:
        # J_next = J + (delta_f-J*delta_x)*delta_x'/(delta_x'delta_x).
        step=alpha.*gradient
        delta_f=ftrial.-f
        correction=delta_f.-J*step
        denominator=dot(step,step)
        denominator>0 || error("vanishing nonlinear step")
        @inbounds for k in 1:n
            coefficient=step[k]/denominator
            for j in 1:n
                J[j,k]+=correction[j]*coefficient
            end
        end
        copyto!(x,xt)
        copyto!(f,ftrial)
        iterations=it
    end
    # Never present a nonconverged path as an impulse response.
    transition_residual!(f,x,input,w)
    throw(ErrorException("transition GE did not converge: iterations=$iterations, "*
                         "max scaled residual=$(maximum(abs,f)), evals=$(w.function_evals)"))
end
