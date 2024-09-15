using HarmonicBalance
using Test

# WARNING: Method definition (::Type{Symbolics.Num})(Base.Complex{Symbolics.Num}) in module HarmonicBalance at E:\HarmonicBalance.jl\src\Symbolics_customised.jl:165 overwritten in module HarmonicBalance on the same line (check for duplicate calls to `include`).
# WARNING: Method definition (::Type{HomotopyContinuation.ModelKit.Variable})(Symbolics.Num) in module HC_wrapper at E:\HarmonicBalance.jl\src\modules\HC_wrapper\homotopy_interface.jl:2 overwritten in module HC_wrapper on the same line (check for duplicate calls to `include`).

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
    using JET
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

@testset "Extentions" begin
    @testset "Time evolution extention" begin
        include("time_evolution.jl")
        include("hysteresis_sweep.jl")
    end
    @testset "ModelingToolkit extention" begin
        include("ModelingToolkitExt.jl")
    end
    @testset "SteadyState Extention" begin
        include("SteadyStateDiffEqExt.jl")
    end
end

@testset "Doctests" begin
    using Documenter
    Documenter.doctest(HarmonicBalance)
end
