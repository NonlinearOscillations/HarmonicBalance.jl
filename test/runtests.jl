using HarmonicBalance
using Test

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)

@testset "Package health" begin
    using ExplicitImports, Aqua
    @test check_no_stale_explicit_imports(HarmonicBalance) == nothing
    @test check_all_explicit_imports_via_owners(HarmonicBalance) == nothing
    Aqua.test_ambiguities(HarmonicBalance)
    # Aqua.test_all(HarmonicBalance)
end

@testset "Symbolics customised" begin
    include("powers.jl")
    include("harmonics.jl")
    include("fourier.jl")
end

@testset "IO" begin
    include("load.jl")
end

@testset "Computing steady states" begin
    include("parametron.jl")
    include("krylov.jl")
end

@testset "Processing solutions" begin
    include("transform_solutions.jl")
end

@testset "Plotting" begin
    include("plotting.jl")
end

@testset "Linear response" begin
    include("linear_response.jl")
end

@testset "Limit cycle" begin
    include("limit_cycle.jl")
end

@testset "Time evolution extention" begin
    include("time_evolution.jl")
    include("hysteresis_sweep.jl")
end

@testset "Extentions" begin
    include("ModelingToolkitExt.jl")
    include("SteadyStateDiffEqExt.jl")
end
