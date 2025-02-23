# # Steady state sweeps

using HarmonicBalance, SteadyStateDiffEq, ModelingToolkit
using BenchmarkTools, Plots, StaticArrays, OrdinaryDiffEq, LinearAlgebra
using HarmonicBalance: OrderedDict

@variables α ω ω0 F γ η t x(t);

diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) + η * d(x, t) * x^2 ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω)
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = (ω0 => 1.0, γ => 1e-2, F => 0.02, α => 1.0, η => 0.3)
ω_span = (0.8, 1.3);
ω_range = range(ω_span..., 201);
varied = ω => ω_range

result_HB = get_steady_states(harmonic_eq, varied, fixed)
plot(result_HB, "sqrt(u1^2+v1^2)")

# ## Steady state sweep using `SteadyStateDiffEq.jl`
param = OrderedDict(merge(Dict(fixed), Dict(ω => 1.1)))
x0 = [0, 0.0];
prob_ss = SteadyStateProblem(
    harmonic_eq, x0, param; in_place=false, u0_constructor=x -> SVector(x...)
)

varied = 6 => ω_range
result_ss = steady_state_sweep(
    prob_ss, DynamicSS(Rodas5()); varied, abstol=1e-5, reltol=1e-5
)

plot(result_HB, "sqrt(u1^2+v1^2)")
plot!(ω_range, norm.(result_ss); c=:gray, ls=:dash)

# ## Adiabatic evolution

timespan = (0.0, 50_000)
sweep = AdiabaticSweep(ω => (0.8, 1.3), timespan) # linearly interpolate between two values at two times
ode_problem = ODEProblem(harmonic_eq, fixed; u0=[0.01; 0.0], timespan, sweep)
time_soln = solve(ode_problem, Tsit5(); saveat=250)

plot(result_HB, "sqrt(u1^2+v1^2)")
plot(time_soln.t, norm.(time_soln.u))

# ## using follow_branch

followed_branch, Ys = follow_branch(1, result_HB; y="√(u1^2+v1^2)")
Y_followed_gr =
    real.([Ys[param_idx][branch] for (param_idx, branch) in enumerate(followed_branch)]);

plot(result_HB, "sqrt(u1^2+v1^2)")
plot!(ω_range, Y_followed_gr; c=:gray, ls=:dash)

# ## comparison

@btime result_ss = steady_state_sweep(
    prob_ss, DynamicSS(Rodas5()); varied, abstol=1e-5, reltol=1e-5
)

@btime time_soln = solve(ode_problem, Tsit5(); saveat=250)

@btime begin
    followed_branch, Ys = follow_branch(1, result_HB; y="√(u1^2+v1^2)")
    Y_followed_gr =
        real.([Ys[param_idx][branch] for (param_idx, branch) in enumerate(followed_branch)])
end

# Plotting them together gives:
plot(ω_range, norm.(result_ss))
plot!(ω_range, norm.(time_soln.u))
plot!(ω_range, Y_followed_gr)
