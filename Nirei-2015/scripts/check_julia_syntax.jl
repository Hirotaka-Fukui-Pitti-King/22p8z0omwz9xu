#!/usr/bin/env julia
# Standalone syntax gate: no external packages or execution of model code.
# Meta.parseall can return Expr(:error) / Expr(:incomplete) inside a
# :toplevel AST rather than raising, so inspect the whole parsed tree.
# Usage from Nirei-2015/: julia --project=. scripts/check_julia_syntax.jl

function first_syntax_error(node)
    node isa Expr || return nothing
    (node.head === :error || node.head === :incomplete) && return node
    for part in node.args
        bad = first_syntax_error(part)
        bad === nothing || return bad
    end
    return nothing
end

function check_source(source::String, filename::String)
    tree = Meta.parseall(source; filename=filename)
    bad = first_syntax_error(tree)
    bad === nothing ||
        error("Julia syntax error in $(filename): $(sprint(show, bad))")
    return nothing
end

function main()
    # Regression for the previously missed reserved-word assignment.
    # A parser accepting this source as valid must fail this gate.
    bad_source = "function sentinel()\n    macro = 1\nend\n"
    bad_tree = Meta.parseall(bad_source; filename="syntax_gate_negative_test.jl")
    first_syntax_error(bad_tree) === nothing &&
        error("syntax gate failed to detect reserved-word assignment")

    root = normpath(joinpath(@__DIR__, ".."))
    count = 0
    for folder in ("src", "scripts")
        dir = joinpath(root, folder)
        for file in sort(readdir(dir; join=true))
            endswith(file, ".jl") || continue
            check_source(read(file, String), file)
            println("Parsed and validated: ", relpath(file, root))
            count += 1
        end
    end
    count > 0 || error("No Julia source files found")
    println("Julia syntax validation passed for ", count, " files.")
    return nothing
end

main()
