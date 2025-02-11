using HarmonicBalance
using ModelingToolkit
using Test

@testset "Utilities" begin
    ModelingToolkitExt = Base.get_extension(HarmonicBalance, :ModelingToolkitExt)
    @variables α
    check = ModelingToolkitExt.declare_parameter(α)
    @test ModelingToolkit.PARAMETER ∈ values(check.val.metadata)
end

@testset "DifferentialEquation" begin
    @testset "ODESystem" begin
        @variables α ω ω0 F γ t x(t)
        diff_eq = DifferentialEquation(
            d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
        )

        fixed = (α => 1.0, ω0 => 1.1, F => 0.01, γ => 0.01)
        param = HarmonicBalance.OrderedDict(merge(Dict(fixed), Dict(ω => 1.1)))
        sys = ODESystem(diff_eq)

        for p in string.([α, ω, ω0, F, γ])
            @test p ∈ string.(parameters(sys))
        end

        # can run a second time without error; diff_eq unmutated
        ODESystem(diff_eq)
    end
    @testset "ODEProblem" begin
        @variables α ω ω0 F γ t x(t)
        diff_eq = DifferentialEquation(
            d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
        )

        add_harmonic!(diff_eq, x, ω) #
        harmonic_eq = get_harmonic_equations(diff_eq)

        sys = ODESystem(harmonic_eq)
        fixed = (α => 1.0, ω0 => 1.1, F => 0.01, γ => 0.01)
        param = HarmonicBalance.OrderedDict(merge(Dict(fixed), Dict(ω => 1.1)))

        ODEProblem(diff_eq, [1.0, 0.0], (0, 100), param)
    end
end

@testset "HarmonicEquation" begin
    @variables α ω ω0 F γ t x(t)
    diff_eq = DifferentialEquation(
        d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
    )

    add_harmonic!(diff_eq, x, ω) #
    harmonic_eq = get_harmonic_equations(diff_eq)

    fixed = (α => 1.0, ω0 => 1.1, F => 0.01, γ => 0.01)
    param = HarmonicBalance.OrderedDict(merge(Dict(fixed), Dict(ω => 1.1)))
    @testset "ODESystem" begin
        sys = ODESystem(harmonic_eq)

        for p in string.([α, ω, ω0, F, γ])
            @test p ∈ string.(parameters(sys))
        end
    end

    @testset "ODEProblem" begin
        ODEProblem(harmonic_eq, [1.0, 0.0], (0, 100), param)
    end
end
