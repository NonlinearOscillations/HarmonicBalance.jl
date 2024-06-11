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
    "krylov.jl",
    "linear_response.jl",
    "limit_cycle.jl"
]

files_ext = [
    "ModelingToolkitExt.jl",
    "SteadyStateDiffEqExt.jl",
    "time_evolution.jl",
    "hysteresis_sweep.jl"
]

# for file in files
#     include(file)
#     printstyled(file * ":    OK\n"; color=:green)
# end

if isdefined(Base, :get_extension) && VERSION >= v"1.9.0"
    for file in files_ext
        include(file)
        printstyled(file * ":    OK\n"; color=:green)
    end
end
printstyled("\nALL TESTS PASSED!\n"; color=:green)
