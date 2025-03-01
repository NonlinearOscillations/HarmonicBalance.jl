using HarmonicBalance, Test

@variables ω₀ γ λ α ω t x(t)

natural_equation = d(d(x, t), t) + γ * d(x, t) + (ω₀^2 - λ * cos(2 * ω * t)) * x + α * x^3
diff_eq = DifferentialEquation(natural_equation, x)

add_harmonic!(diff_eq, x, ω);

harmonic_eq = get_harmonic_equations(diff_eq)

fixed = (ω₀ => 1.0, γ => 0.002, α => 1.0)
varied = (ω => range(0.99, 1.01, 100), λ => range(1e-6, 0.03, 100))

result_2D = get_steady_states(harmonic_eq, varied, fixed; verbose=true)

@testset "filter solutions" begin
    result = deepcopy(result_2D)
    @test length(result.solutions[1]) == 5
    filter_result!(result, "stable")
    @test length(result.solutions[1]) == 3
end

@testset "Classifying solutions" begin
    result = deepcopy(result_2D)
    filter_result!(result, "stable")
    classify_solutions!(result, "sqrt(u1^2 + v1^2) >0.0", "not_zero")
    @test sum(any.(get_class(result, "not_zero"))) == 2
end
