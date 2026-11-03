
# No standalone plotting of guessed or unverified paths.
# The stage3 run validates BOTH heterogeneous-agent regimes before
# generating all three SVG figures.
include(joinpath(@__DIR__,"run_sudden_stop.jl"))
