using HarmonicBalance

@variables Ω γ λ F x θ η α ω0 ω t T ψ
@variables x(t) y(t)

natural_equation = [
    d(d(x, t), t) + γ * d(x, t) + Ω^2 * x + α * x^3 ~ F * cos(ω * t),
    d(d(y, t), t) + γ * d(y, t) + Ω^2 * y + α * y^3 ~ 0,
]
dEOM = DifferentialEquation(natural_equation, [x, y])

add_harmonic!(dEOM, x, ω)
add_harmonic!(dEOM, y, ω)
harmonic_eq = get_harmonic_equations(dEOM);

fixed = (Ω => 1.0, γ => 1e-2, F => 1e-3, α => 1.0)
varied = ω => range(0.9, 1.1, 10)
res = get_steady_states(harmonic_eq, varied, fixed; show_progress=false);

transform_solutions(res, "u1^2+v1^2")
transform_solutions(res, "√(u1^2+v1^2)"; realify=true)

@testset "to_lab_frame" begin
    using HarmonicBalance: to_lab_frame
    @variables z(t)
    times = 0:1:10
    @test to_lab_frame(res, x, times; index=1, branch=1) != zeros(length(times))
    @test all(isapprox.(to_lab_frame(res, y, times; index=1, branch=1), 0.0, atol=1e-10))
    @test all(to_lab_frame(res, z, times; index=1, branch=1) .≈ zeros(length(times)))
    @test to_lab_frame(res, d(x, t), times; index=1, branch=1) != zeros(length(times))
end

using HarmonicBalance
