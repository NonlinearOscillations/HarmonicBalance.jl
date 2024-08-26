using HarmonicBalance
using Symbolics, OrdinaryDiffEq
using Test

@variables Ω, γ, λ, F, x, θ, η, α, ω0, ω, t, T, ψ
@variables x(t)

natural_equation =
    d(d(x, t), t) +
    γ * d(x, t) +
    Ω^2 * (1 - λ * cos(2 * ω * t + ψ)) * x +
    α * x^3 +
    η * d(x, t) * x^2
forces = F * cos(ω * t + θ)
dEOM = DifferentialEquation(natural_equation + forces, x)
add_harmonic!(dEOM, x, ω)
harmonic_eq = get_harmonic_equations(dEOM; slow_time=T, fast_time=t);
p = HarmonicBalance.Problem(harmonic_eq);

fixed = (Ω => 1.0, γ => 1e-2, λ => 5e-2, F => 1e-3, α => 1.0, η => 0.3, θ => 0, ψ => 0)
varied = ω => range(0.9, 1.1, 10)
res = get_steady_states(p, varied, fixed; show_progress=false, seed=SEED)

tspan = (0.0, 10)
sweep = ParameterSweep(ω => (0.9, 1.1), tspan) # linearly interpolate between two values at two times
ode_problem = ODEProblem(harmonic_eq, fixed; sweep=sweep, x0=[0.01; 0.0], timespan=tspan)
time_soln = solve(ode_problem, Tsit5(); saveat=1);

transform_solutions(time_soln, "sqrt(u1^2+v1^2)", harmonic_eq)
