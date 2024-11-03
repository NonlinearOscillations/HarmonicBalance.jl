using HarmonicBalance, OrdinaryDiffEqTsit5

@variables α ω ω0 F γ η t x(t); # declare constant variables and a function x(t)

diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0 * x + α * x^3 + γ * d(x, t) + η * x^2 * d(x, t) ~ F * cos(ω * t), x
) # define ODE

add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)
harmonic_eq = get_harmonic_equations(diff_eq) # implement ansatz to get harmonic equations

fixed = (α => 1, ω0 => 1.0, γ => 0.005, F => 0.005, η => 0.2)   # fixed parameters
varied = ω => range(0.95, 1.1, 10)           # range of parameter values
method = HarmonicBalance.WarmUp(; seed=SEED)
result = get_steady_states(harmonic_eq, method, varied, fixed; show_progress=false)

followed_branches, _ = follow_branch(1, result)
followed_branches, _ = follow_branch(1, result; y="√(u1^2+v1^2)")

@test first(followed_branches) ≠ last(followed_branches)

plot_1D_solutions_branch(1, result; x="ω", y="u1^2+v1^2")

plot_linear_response(result, x, followed_branches; Ω_range=range(0.9, 1.1, 10), show=false);

# Check if an empty brspectrum can be plotted
followed_branches[6] = 3
plot_linear_response(
    result, x, followed_branches; Ω_range=range(0.9, 1.1, 10), show=false, force=true
);
