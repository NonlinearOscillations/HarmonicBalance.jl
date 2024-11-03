using HarmonicBalance, Plots
default(; show=false)
#using Test #

@variables γ λ x η α ω0 ω
@variables t x(t)

natural_equation =
    d(d(x, t), t) +
    γ * d(x, t) +
    ω0^2 * (1 - λ * cos(2 * ω * t)) * x +
    α * x^3 +
    η * d(x, t) * x^2
dEOM = DifferentialEquation(natural_equation, x)
add_harmonic!(dEOM, x, ω)
harmonic_eq = get_harmonic_equations(dEOM);

fixed = (ω0 => 1.0, γ => 1e-2, λ => 5e-2, α => 1.0, η => 0.3)
varied = ω => range(0.9, 1.1, 10)
res = get_steady_states(harmonic_eq, varied, fixed; show_progress=false)

# plot 1D result
plot(res; x="ω", y="u1");
plot(res; y="√(u1^2+v1^2)");
plot_spaghetti(res; x="v1", y="u1", z="ω");
plot_phase_diagram(res);

fixed = (ω0 => 1.0, γ => 1e-2, α => 1.0, η => 0.3)
varied = (ω => range(0.9, 1.1, 5), λ => range(0.01, 0.05, 5))
res = get_steady_states(harmonic_eq, varied, fixed; show_progress=false)

# plot 2D result
plot_phase_diagram(res);
plot(res, "√(u1^2+v1^2)"; branch=1);
plot(res; y="√(u1^2+v1^2)", cut=λ => 0.03);
@test_throws ArgumentError plot(res; y="√(u1^2+v1^2)", cut=λ => 0.1);
