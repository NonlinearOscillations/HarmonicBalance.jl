using Test
using HarmonicBalance
@variables α ω ω0 F γ η t x(t); # declare constant variables and a function x(t)

diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0 * x + α * x^3 + γ * d(x, t) + η * x^2 * d(x, t) ~ F * cos(ω * t), x
) # define ODE
add_harmonic!(diff_eq, x, ω)
harmonic_eq = get_harmonic_equations(diff_eq)

@testset "default variable names" begin
    vars = get_variables(harmonic_eq)
    @test HarmonicBalance.var_name.(vars) == ["u1", "v1"]
end
