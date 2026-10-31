#!/usr/bin/env julia
# Standalone syntax gate: no external packages and no execution of model code.
# Usage from Nirei-2015/: julia --project=. scripts/check_julia_syntax.jl

function main()
    root = normpath(joinpath(@__DIR__, ".."))
    count = 0
    for folder in ("src", "scripts")
        dir = joinpath(root, folder)
        for file in sort(readdir(dir; join=true))
            endswith(file, ".jl") || continue
            Meta.parseall(read(file, String); filename=file)
            println("Parsed: ", relpath(file, root))
            count += 1
        end
    end
    count > 0 || error("No Julia source files found")
    println("Julia syntax validation passed for ", count, " files.")
    return nothing
end

main()
