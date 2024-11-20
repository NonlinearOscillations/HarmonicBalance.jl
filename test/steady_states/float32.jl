using HarmonicBalance
using Test

@testset "Float32" begin
    @variables ω₀ γ λ α ω t x(t)
    natural_equation =
        d(d(x, t), t) + γ * d(x, t) + (ω₀^2 - λ * cos(2 * ω * t)) * x + α * x^3
    diff_eq = DifferentialEquation(natural_equation, x)

    add_harmonic!(diff_eq, x, ω)

    harmonic_eq = get_harmonic_equations(diff_eq)

    fixed = (ω₀ => 1.0f0, γ => 0.005f0, λ => 0.05f0, α => 1.0f0)
    varied = ω => range(0.9f0, 1.1f0, 100)

    method = WarmUp{ComplexF32}()
    result = get_steady_states(harmonic_eq, method, varied, fixed)
    @test typeof(result).parameters[1] == ComplexF32 broken = true
end
