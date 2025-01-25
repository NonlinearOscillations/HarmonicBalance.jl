|
@testset "Concretely typed" begin
    CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
    julia_version = VERSION >= v"1.11.0-DEV.0" # fails on 1.11
    if !CI && !julia_version
        using HarmonicBalance

        using CheckConcreteStructs

        all_concrete(HarmonicBalance.HarmonicVariable)
        all_concrete(HarmonicBalance.WarmUp)
        all_concrete(HarmonicBalance.TotalDegree)
        all_concrete(HarmonicBalance.Polyhedral)
        all_concrete(HarmonicBalance.Result)
        all_concrete(HarmonicBalance.HomotopyContinuationProblem)
        all_concrete(HarmonicBalance.HarmonicEquation)
        all_concrete(HarmonicBalance.DifferentialEquation)
        all_concrete(HarmonicBalance.HarmonicVariable)
        all_concrete(HarmonicBalance.AdiabaticSweep)

        all_concrete(HarmonicBalance.LinearResponse.Lorentzian)
        all_concrete(HarmonicBalance.LinearResponse.ResponseMatrix)
        all_concrete(HarmonicBalance.LinearResponse.JacobianSpectrum)
    end
end

@testset "Code linting" begin
    using JET
    JET.test_package(HarmonicBalance; target_defined_modules=true)
end

@testset "Code quality" begin
    using ExplicitImports, Aqua
    using ModelingToolkit, OrdinaryDiffEqTsit5, SteadyStateDiffEq
    ignore_deps = [:Random, :LinearAlgebra, :Printf, :Test, :Pkg]
    TimeEvolution = Base.get_extension(HarmonicBalance, :TimeEvolution)
    ModelingToolkitExt = Base.get_extension(HarmonicBalance, :ModelingToolkitExt)
    SteadyStateDiffEqExt = Base.get_extension(HarmonicBalance, :SteadyStateDiffEqExt)
    @test check_no_stale_explicit_imports(HarmonicBalance) == nothing
    @test check_all_explicit_imports_via_owners(HarmonicBalance) == nothing
    Aqua.test_ambiguities([HarmonicBalance])
    Aqua.test_all(
        HarmonicBalance;
        deps_compat=(
            ignore=ignore_deps,
            check_extras=(ignore=ignore_deps,),
            check_weakdeps=(ignore=ignore_deps,),
        ),
        ambiguities=false,
    )
    for mod in [TimeEvolution, ModelingToolkitExt, SteadyStateDiffEqExt]
        @test check_no_stale_explicit_imports(mod) == nothing
        @test check_all_explicit_imports_via_owners(mod) == nothing
        # Aqua.test_ambiguities(mod)
        Aqua.test_all(
            mod;
            deps_compat=false,
            ambiguities=false,
            piracies=false,
            stale_deps=false,
            project_extras=false,
            persistent_tasks=false,
        )
    end
end
