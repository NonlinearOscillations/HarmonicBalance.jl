using HarmonicBalance
using BenchmarkTools

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)
@variables t x(t) y(t) z(t);
@variables ω0 ω γ λ α η J F;

equations = [
    d(d(x, t), t) + ω0^2 * (1 - λ * cos(2 * ω * t)) * x + γ * d(x, t) + α * x^3 - J * y -
    J * z,
    d(d(y, t), t) + ω0^2 * (1 - λ * cos(2 * ω * t)) * y + γ * d(y, t) + α * y^3 - J * x -
    J * z,
    d(d(z, t), t) + ω0^2 * (1 - λ * cos(2 * ω * t)) * z + γ * d(z, t) + α * z^3 - J * x -
    J * y,
]

system = DifferentialEquation(equations, [x, y, z])

add_harmonic!(system, x, ω) # x will rotate at ω
add_harmonic!(system, y, ω) # y will rotate at ω
add_harmonic!(system, z, ω) # z will rotate at ω

harmonics = get_harmonic_equations(system);

prob = HarmonicBalance.Problem(harmonics)

res = 10
fixed = Dict(ω0 => 1.0, γ => 0.005, α => 1.0, η => 0.0, J => 0.005)
varied = (ω => range(0.985, 1.015, res), λ => range(1e-6, 0.04, res))

@btime res = get_steady_states(prob, WarmUp(), varied, fixed; show_progress=false) # 399.727 ms (561289 allocations: 63.11 MiB)
@btime res = get_steady_states(
    prob, WarmUp(; compile=true), varied, fixed; show_progress=false
)
@btime res = get_steady_states(
    prob, Polyhedral(; only_non_zero=true), varied, fixed; show_progress=false
)
@btime res = get_steady_states(
    prob, Polyhedral(; only_non_zero=false), varied, fixed; show_progress=false
) # 382.672 ms (930906 allocations: 51.16 MiB)
@btime res = get_steady_states(
    prob, TotalDegree(; compile=true), varied, fixed; show_progress=false
) # 378.314 ms (904817 allocations: 49.63 MiB)
@btime res = get_steady_states(prob, TotalDegree(), varied, fixed; show_progress=false)# 379.858 ms (925130 allocations: 50.86 MiB)
plot(res; y="u1")
