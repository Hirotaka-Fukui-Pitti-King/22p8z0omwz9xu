
"""
The transition household block follows author backsolve_egm:
given the terminal steady-state consumption policy, move t=T-1,...,1
backwards. Euler uses next-period gross return 1+r[t+1], while the
current budget uses 1+r[t]. Fine-grid policies are evaluated DIRECTLY
on the EGM endogenous current-asset nodes.
"""
function transition_household_backward!(w::TransitionWorkspace{T},
                                        input::TransitionInputs{T}) where {T}
    p=input.params
    coarse=input.grids.coarse
    fine=input.grids.fine
    chain=input.chain
    nc=length(coarse)
    nf=length(fine)
    ns=length(chain.levels)
    n=length(input.B)
    copyto!(view(w.coarse_consumption,:,:,n),
            input.terminal.household.consumption_coarse)
    copyto!(view(w.fine_consumption,:,:,n),
            input.terminal.household.consumption_fine)
    copyto!(view(w.fine_savings,:,:,n),
            input.terminal.household.savings_fine)
    @inbounds for t in (n-1):-1:1
        now=view(w.coarse_consumption,:,:,t)
        future=view(w.coarse_consumption,:,:,t+1)
        cnow=view(w.fine_consumption,:,:,t)
        anow=view(w.fine_savings,:,:,t)
        gross_old=1+w.rate[t]
        gross_next=1+w.rate[t+1]
        m=consumption_cost(w.terms[t],p)
        m>0 || throw(DomainError(m,"nonpositive expenditure conversion"))
        for j in 1:ns
            # Gross labor earnings equal author markup*(1-alpha)*Y*ex,
            # since wage*L = markup*(1-alpha)*Y.
            income=w.wage[t]*w.L[t]*chain.levels[j]+w.profits[t]
            for i in 1:nc
                emu=zero(T)
                for jp in 1:ns
                    emu+=chain.transition[j,jp]*future[i,jp]^(-p.gamma)
                end
                c_endo=(input.initial.beta*gross_next*emu)^(-inv(p.gamma))
                w.endogenous_consumption[i]=c_endo
                w.endogenous_asset[i]=(m*c_endo+coarse[i]-income)/gross_old
            end
            all(diff(w.endogenous_asset) .> 0) ||
                throw(ArgumentError("non-monotone backward EGM endogenous assets at t=$t,s=$j"))
            pchip_slopes!(w.interpolation_derivative,
                          w.endogenous_asset,w.endogenous_consumption)
            cutoff=w.endogenous_asset[1]
            for i in 1:nc
                a=coarse[i]
                c= if a<=cutoff
                    (gross_old*a+income-input.grids.borrowing_limit)/m
                else
                    pchip_evaluate(w.endogenous_asset,w.endogenous_consumption,
                                   w.interpolation_derivative,a)
                end
                isfinite(c) && c>0 ||
                    throw(DomainError(c,"nonpositive transition coarse consumption at t=$t"))
                now[i,j]=c
            end
            for i in 1:nf
                a=fine[i]
                if a<=cutoff
                    saving=input.grids.borrowing_limit
                    cons=(gross_old*a+income-saving)/m
                else
                    cons=pchip_evaluate(w.endogenous_asset,w.endogenous_consumption,
                                        w.interpolation_derivative,a)
                    saving=gross_old*a+income-m*cons
                end
                isfinite(cons) && cons>0 ||
                    throw(DomainError(cons,"nonpositive transition fine consumption at t=$t"))
                isfinite(saving) && saving>=input.grids.borrowing_limit ||
                    throw(DomainError(saving,"transition borrowing violation at t=$t"))
                # The author backsolve_egm checks the borrowing
                # restriction here, not the upper distribution-grid bound.
                # Endogenous policy may leave the grid at far-right states
                # carrying ZERO probability. Retain the true saving value;
                # distribution_step! rejects out-of-grid choices only
                # when a state actually has positive probability mass.
                cnow[i,j]=cons
                anow[i,j]=saving
            end
        end
    end
    return w
end
