using HarmonicBalance
using Symbolics, Plots
default(show=false)
using Test

@variables γ, λ, x, η, α, ω0, ω
@variables t x(t)

natural_equation =  d(d(x,t),t) + γ*d(x,t) + ω0^2*(1-λ*cos(2*ω*t))*x + α*x^3 + η*d(x,t)*x^2
dEOM = DifferentialEquation(natural_equation, x)
add_harmonic!(dEOM, x, ω)
harmonic_eq = get_harmonic_equations(dEOM);

fixed = (ω0 => 1.0, γ => 1e-2, λ => 5e-2, α => 1.0, η => 0.3)
varied = ω => range(0.9, 1.1, 100)
res = get_steady_states(harmonic_eq, varied, fixed)

# plot 1D result
plot(res, x="ω", y="u1");

fixed = (ω0 => 1.0, γ => 1e-2, α => 1.0, η => 0.3)
varied = (ω => LinRange(0.9, 1.1, 10), λ => LinRange(0.01, 0.05, 10))
res = get_steady_states(harmonic_eq, varied, fixed)

# plot 2D result
plot_phase_diagram(res);


