# Harmonic oscillator: comparison of KB and HB methods

using HarmonicBalance

@variables ω ω0 F γ t x(t) # declare constant variables and a function x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)

krylov_eq1 = get_krylov_equations(diff_eq; order=1)
krylov_eq2 = get_krylov_equations(diff_eq; order=2)
harmonic_eq = get_harmonic_equations(diff_eq)
harmonic_eq = rearrange_standard(harmonic_eq)

varied = (ω => range(0.1, 1.9, 200)) # range of parameter values
fixed = (ω0 => 1.0, γ => 0.05, F => 0.1) # fixed parameters
result_krylov1 = get_steady_states(krylov_eq1, varied, fixed)
result_krylov2 = get_steady_states(krylov_eq2, varied, fixed)
result_harmonic = get_steady_states(harmonic_eq, varied, fixed)

plot(result_krylov1, y="u1^2+v1^2", label="Krylov 1")
plot(result_krylov2, y="u1^2+v1^2", label="Krylov 2")
plot(result_harmonic, y="u1^2+v1^2", label="Harmonic")

plot(result_krylov1, y="u1", label="Krylov 1", legend=:best)
plot(result_krylov2, y="u1", label="Krylov 2", legend=:best)
plot(result_harmonic, y="u1", label="Harmonic", legend=:best)

plot(result_krylov1, y="v1", label="Krylov 1", legend=:best)
plot(result_krylov2, y="v1", label="Krylov 2", legend=:best)
plot(result_harmonic, y="v1", label="Harmonic", legend=:best)

plot_eigenvalues(result_krylov1, branch=1)
plot_eigenvalues(result_krylov2, branch=1)
plot_eigenvalues(result_harmonic, branch=1)

plot_linear_response(result_krylov1, x; Ω_range=range(0.1, 1.9, 200), branch=1)
plot_linear_response(result_krylov2, x; Ω_range=range(0.1, 1.9, 200), branch=1)
plot_linear_response(result_harmonic, x; Ω_range=range(0.1, 1.9, 200), branch=1)
