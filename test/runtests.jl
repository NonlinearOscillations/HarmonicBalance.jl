using Pkg
current_path = @__DIR__
Pkg.activate(current_path * "/../.");
using HarmonicBalance
using Test

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)

files = [
    "powers.jl",
    "harmonics.jl",
    "fourier.jl",
    "load.jl",
    "parametron.jl",
    "transform_solutions.jl",
    "plotting.jl",
    "time_evolution.jl",
    "krylov.jl",
    "hysteresis_sweep.jl",
    "linear_response.jl",
    "limit_cycle.jl",
    "ModelingToolkitExt.jl"
]

for file in files
    include(file)
    printstyled(file * ":    OK\n"; color=:green)
end

printstyled("\nALL TESTS PASSED!\n"; color=:green)
