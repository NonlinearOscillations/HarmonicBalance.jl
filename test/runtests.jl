using HarmonicBalance
using Test

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)

@testset "Code quality" begin
    using ExplicitImports, Aqua
    ignore_deps = [:Random, :LinearAlgebra, :Printf, :Test, :Pkg]

    @test check_no_stale_explicit_imports(HarmonicBalance) == nothing
    @test check_all_explicit_imports_via_owners(HarmonicBalance) == nothing
    Aqua.test_ambiguities(HarmonicBalance)
    Aqua.test_all(
        HarmonicBalance;
        deps_compat=(
            ignore=ignore_deps,
            check_extras=(ignore=ignore_deps,),
            check_weakdeps=(ignore=ignore_deps,),
        ),
        piracies=(treat_as_own=[HarmonicBalance.Num],),
        ambiguities=false,
    )
end

@testset "Code linting" begin
    JET.test_package(HarmonicBalance; target_defined_modules=true)
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

@testset "Doctests" begin
    Documenter.doctest(HiddenMarkovModels)
end
