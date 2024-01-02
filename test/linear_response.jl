using HarmonicBalance

import HarmonicBalance.LinearResponse.plot_linear_response
import HarmonicBalance.LinearResponse.plot_rotframe_jacobian_response
import HarmonicBalance.LinearResponse.plot_eigenvalues

@variables α, ω, ω0, F, γ, t, x(t);

diff_eq = DifferentialEquation(d(x,t,2) + ω0*x + α*x^3 + γ*d(x,t) ~ F*cos(ω*t), x)
add_harmonic!(diff_eq, x, ω)
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = (α => 1, ω0 => 1.0, γ => 1e-2, F => 1e-6)
varied = ω => range(0.9, 1.1, 10)

result = get_steady_states(harmonic_eq, varied, fixed, show_progress=false, seed=SEED)

plot_linear_response(result, x, branch=1, Ω_range=range(0.9,1.1,10), order=1, logscale=true)

plot_rotframe_jacobian_response(result, Ω_range=range(0.01,1.1,10), branch=1, logscale=true)
