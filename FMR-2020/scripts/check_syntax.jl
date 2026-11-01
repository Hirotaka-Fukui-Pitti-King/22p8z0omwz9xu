# Run from the project with: julia --project=. scripts/check_syntax.jl
root_dir = normpath(joinpath(@__DIR__,".."))
count = 0
for (dir, _, names) in walkdir(root_dir)
    for name in names
        endswith(name,".jl") || continue
        path = joinpath(dir,name)
        # Meta.parseall may return an AST with nested Expr(:error) or
        # Expr(:incomplete) instead of throwing a ParseError. Explicitly
        # reject those expressions; printing "PARSE OK" is insufficient.
        expr = Meta.parseall(read(path,String);filename=path)
        function check_parse_tree(node)
            if node isa Expr
                (node.head === :error || node.head === :incomplete) &&
                    error("Julia parse failed for $(path): $(sprint(show,node))")
                for item in node.args
                    check_parse_tree(item)
                end
            end
            return nothing
        end
        check_parse_tree(expr)
        println("PARSE OK: ",relpath(path,root_dir))
        global count += 1
    end
end
count > 0 || error("No Julia files found")
println("Parsed ",count," Julia files (parse only; does not execute).")
