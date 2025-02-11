# Steady state sweeps

using HarmonicBalance, SteadyStateDiffEq, ModelingToolkit
using BenchmarkTools, Plots, StaticArrays, OrdinaryDiffEq, LinearAlgebra
using HarmonicBalance: OrderedDict

@variables α ω ω0 F γ t x(t);

diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω)
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = (ω0 => 1.0, γ => 0.005, F => 0.0052, α => 1.0)
ω_span = (0.8, 1.5);
ω_range = range(ω_span..., 200);
varied = ω => ω_range

result_HB = get_steady_states(harmonic_eq, varied, fixed)
plot(result_HB, "sqrt(u1^2+v1^2)")

## Steady state sweep using `SteadyStateDiffEq.jl`
param = OrderedDict(merge(Dict(fixed), Dict(ω => 1.1)))
x0 = [1.0, 0.0];
prob_ss = SteadyStateProblem(
    harmonic_eq, x0, param; in_place=false, u0_constructor=x -> SVector(x...)
)

sweep = 5 => ω_range
result_ss = steady_state_sweep(
    prob_ss, DynamicSS(Rodas5()); varied=sweep, abstol=1e-5, reltol=1e-5
)

plot(result_HB, "sqrt(u1^2+v1^2)")
plot!(ω_range, norm.(result_ss))

## Adiabatic evolution

timespan = (0.0, 10_000)
sweep = AdiabaticSweep(ω => (0.8, 1.3), tspan) # linearly interpolate between two values at two times
ode_problem = ODEProblem(harmonic_eq, fixed; u0=[0.01; 0.0], timespan, sweep)
time_soln = solve(ode_problem, Tsit5(); saveat=100, dt=1e-5)

plot(result_HB, "sqrt(u1^2+v1^2)")
plot(time_soln.t, norm.(time_soln.u))

## using follow_branch

followed_branch, Ys = follow_branch(1, result_HB; y="√(u1^2+v1^2)")
Y_followed_gr =
    real.([Ys[param_idx][branch] for (param_idx, branch) in enumerate(followed_branch)]);

plot(result_HB, "sqrt(u1^2+v1^2)")
plot!(ω_range, Y_followed_gr)
