using HarmonicBalance

using Random
const SEED = 0xd8e5d8df
Random.seed!(SEED)

@variables Ω γ λ F x θ η α ω0 ω t T ψ
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

@time harmonic_eq = get_harmonic_equations(dEOM; slow_time=T, fast_time=t);
@time prob = HarmonicBalance.Problem(harmonic_eq)

fixed = (Ω => 1.0, γ => 1e-2, λ => 5e-2, F => 0, α => 1.0, η => 0.3, θ => 0, ψ => 0)
varied = ω => range(0.9, 1.1, 20)
@time res = get_steady_states(prob, varied, fixed; show_progress=false)
