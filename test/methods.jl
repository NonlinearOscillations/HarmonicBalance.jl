using HarmonicBalance
using Test

@testset "WarmUp" begin
    @variables α ω ω0 λ γ t x(t)
    diff_eq = DifferentialEquation(
        d(d(x, t), t) + γ * d(x, t) + ω0^2 * (1 - λ * cos(2 * ω * t)) * x + α * x^3, x
    )

    add_harmonic!(diff_eq, x, ω) #
    harmonic_eq = get_harmonic_equations(diff_eq)

    fixed = HarmonicBalance.OrderedDict((α => 1.0, ω0 => 1.1, λ => 0.01, γ => 0.01))
    varied = HarmonicBalance.OrderedDict((ω => range(0.9, 1.1, 20),))
    prob = HarmonicBalance.Problem(harmonic_eq)

    unique_fixed, input_array = HarmonicBalance._prepare_input_params(prob, varied, fixed)
    @test length.(input_array) == fill(5, 20)
end

@testset "Polyhedral" begin
    @variables α ω ω0 λ γ t x(t)
    diff_eq = DifferentialEquation(
        d(d(x, t), t) + γ * d(x, t) + ω0^2 * (1 - λ * cos(2 * ω * t)) * x + α * x^3, x
    )

    add_harmonic!(diff_eq, x, ω) #
    harmonic_eq = get_harmonic_equations(diff_eq)

    fixed = HarmonicBalance.OrderedDict((α => 1.0, ω0 => 1.1, λ => 0.001, γ => 0.01))
    varied = HarmonicBalance.OrderedDict((ω => range(0.9, 1.1, 10),))
    result = get_steady_states(
        harmonic_eq, Polyhedral(; only_non_zero=true), varied, fixed; show_progress=false
    )
    @test sum(reduce(vcat, classify_branch(result, "stable"))) == 0
end
