using HarmonicBalance
using ModelingToolkit
using ModelingToolkit: varmap_to_vars
using Test

@variables α ω ω0 F γ t x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω) #
harmonic_eq = get_harmonic_equations(diff_eq)

sys = ODESystem(harmonic_eq)
fixed = (α => 1.0, ω0 => 1.1, F => 0.01, γ => 0.01)
param = ParameterList(merge(Dict(fixed), Dict(ω => 1.1)))

for p in string.([α, ω, ω0, F, γ])
    @test p ∈ string.(parameters(sys))
end

ODEProblem(harmonic_eq, [1.0, 0.0], (0, 100), param)
