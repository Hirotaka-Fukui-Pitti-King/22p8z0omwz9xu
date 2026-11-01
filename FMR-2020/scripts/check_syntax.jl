# Run from the project with: julia --project=. scripts/check_syntax.jl
root_dir = normpath(joinpath(@__DIR__,".."))
count = 0
for (dir, _, names) in walkdir(root_dir)
    for name in names
        endswith(name,".jl") || continue
        path = joinpath(dir,name)
        Meta.parseall(read(path,String))
        println("PARSE OK: ",relpath(path,root_dir))
        global count += 1
    end
end
count > 0 || error("No Julia files found")
println("Parsed ",count," Julia files (parse only; does not execute).")
