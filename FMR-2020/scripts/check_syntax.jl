# Run from the project with: julia --project=. scripts/check_syntax.jl
root_dir = normpath(joinpath(@__DIR__,".."))
count = 0
# Parsing does not catch two adjacent Julia docstrings: Julia attaches the
# first to the second string and package loading then fails with "cannot
# document the following expression". Reject that exact pattern EARLY.
triple_quote = "\"\"\""
adjacent_docstrings = Regex("(?s)" * triple_quote * ".*?" *
                            triple_quote * "\\s*" * triple_quote)
example_bad = triple_quote * "first" * triple_quote * "\n" *
              triple_quote * "second" * triple_quote
occursin(adjacent_docstrings, example_bad) ||
    error("Adjacency guard self-test failed")
for (dir, _, names) in walkdir(root_dir)
    for name in names
        endswith(name,".jl") || continue
        path = joinpath(dir,name)
        # Meta.parseall may return an AST with nested Expr(:error) or
        # Expr(:incomplete) instead of throwing a ParseError. Explicitly
        # reject those expressions; printing "PARSE OK" is insufficient.
        source = read(path,String)
        occursin(adjacent_docstrings,source) &&
            error("Consecutive Julia docstrings in $(relpath(path,root_dir)); " *
                  "replace an unattached documentation string with # comments")
        expr = Meta.parseall(source;filename=path)
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
println("Parsed ",count," Julia files (parse only).")

# A valid Meta.parseall result does not prove that Julia can evaluate the
# module: malformed docstring attachment is one example. Load the actual
# package after parsing every Julia file to catch such top-level errors.
using FMR2020
println("LOAD OK: FMR2020 (module loaded; economic tests not executed)")
