# Export the ACTUAL Julia calibrated steady-state prices and terminal
# 161-node policy, together with the original 193/769 full-grid inputs.
# This is an input-only fixture; no 399-variable solve is attempted.
using FMR2020

function write_numeric_csv(path::AbstractString, data::AbstractArray)
    mat=ndims(data)==1 ? reshape(data,:,1) : data
    ndims(mat)==2 || throw(DimensionMismatch("CSV fixture must be 1D/2D"))
    open(path,"w") do io
        for i in axes(mat,1)
            println(io,join((repr(mat[i,j]) for j in axes(mat,2)),","))
        end
    end
    return path
end

function main()
    length(ARGS)==1 || error("Pass calibrated fixture directory")
    dir=ARGS[1]
    mkpath(dir)
    states=solve_three_steady_states(default_parameters();verbose=false)
    st=states.calibrated
    grid=states.grids
    ch=states.chain
    values=Float64[st.beta,st.r,st.terms,st.wage,st.L,st.profits,
                   states.parameters.gamma,states.parameters.theta,
                   one(st.beta)-states.parameters.chi]
    write_numeric_csv(joinpath(dir,"prices.csv"),values)
    write_numeric_csv(joinpath(dir,"coarse_grid.csv"),grid.coarse_full)
    write_numeric_csv(joinpath(dir,"fine_grid.csv"),grid.fine_full)
    write_numeric_csv(joinpath(dir,"productivity.csv"),ch.levels)
    write_numeric_csv(joinpath(dir,"transition.csv"),ch.transition)
    write_numeric_csv(joinpath(dir,"julia_active_coarse_policy.csv"),
                      st.household.consumption_coarse)
    write_numeric_csv(joinpath(dir,"julia_active_fine_policy.csv"),
                      st.household.consumption_fine)
    println("REAL CALIBRATED INPUT: beta=",st.beta," r=",st.r,
            " L=",st.L," coarse_full=",length(grid.coarse_full),
            " coarse_active=",length(grid.coarse),
            " fine_full=",length(grid.fine_full),
            " fine_active=",length(grid.fine))
    println("CALIBRATED INPUT EXPORTED (not a source GE solution): ",dir)
end

main()
