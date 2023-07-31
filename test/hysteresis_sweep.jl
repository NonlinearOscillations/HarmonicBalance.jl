using HarmonicBalance, OrdinaryDiffEq

@variables α ω ω0 F γ η t x(t); # declare constant variables and a function x(t)

diff_eq = DifferentialEquation(d(x,t,2) + ω0*x + α*x^3 + γ*d(x,t) + η*x^2*d(x,t) ~ F*cos(ω*t), x) # define ODE

add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)
harmonic_eq = get_harmonic_equations(diff_eq) # implement ansatz to get harmonic equations

fixed = (α => 1, ω0 => 1.0, γ => 0.005, F => 0.005, η => 0.2)   # fixed parameters
varied = ω => range(0.95, 1.1, 10)           # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed, show_progress=false)

followed_branch, _ = follow_branch(1, result)

@test first(followed_branch) ≠ last(followed_branch)
