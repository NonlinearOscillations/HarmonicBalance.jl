using HarmonicBalance
using Test

using Random
const SEED = 0x8f88209c
Random.seed!(SEED)

@testset "Code quality" begin
    include("code_quality.jl")
end

@testset "API" begin
    include("API.jl")
    include("Problem.jl")
end

@testset "Symbolics customised" begin
    include("symbolics.jl")
end

@testset "IO" begin
    include("load.jl")
end

@testset "Computing steady states" begin
    include("parametron.jl")
    include("krylov.jl")
    include("methods.jl")
end

@testset "Processing solutions" begin
    include("Jacobian.jl")
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

@testset "extensions" begin
    @testset "Time evolution extension" begin
        include("time_evolution.jl")
        include("hysteresis_sweep.jl")
    end
    @testset "ModelingToolkit extension" begin
        include("ModelingToolkitExt.jl")
    end
    @testset "SteadyState extension" begin
        include("SteadyStateDiffEqExt.jl")
    end
end

@testset "Doctests" begin
    using Documenter
    Documenter.doctest(HarmonicBalance)
end
