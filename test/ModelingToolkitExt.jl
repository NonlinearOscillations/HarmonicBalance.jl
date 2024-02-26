using HarmonicBalance
using ModelingToolkit
using Test

@variables α ω ω0 F γ t x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x)
add_harmonic!(diff_eq, x, ω) #
harmonic_eq = get_harmonic_equations(diff_eq)

ODESystem(harmonic_eq)

force = 0.01
omega0 = 1.1
alpha = 1.0
gamma = 0.01
fixed_nonlin = (α => alpha, ω0 => omega0, F => force, γ => gamma)
ω_span = (0.9, 1.5)
ω_range = range(ω_span..., 100)
varied = ω => ω_range

param = ParameterList(merge(Dict(fixed_nonlin), Dict(ω => 1.1)))
varied = 1 => ω_range
x0 = [1.0, 0.0]

ODEProblem(harmonic_eq, x0, (0, 100), param)
