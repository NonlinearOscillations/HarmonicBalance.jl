using HarmonicBalance
using QuantumCumulants

using QuantumCumulants, HarmonicBalance
using LinearAlgebra, Symbolics

@variables α, ωₚ, ω0, t, η, x(t), y(t)

# differential equations
diff_eq = DifferentialEquation([d(x, t, 2) + ω0^2 * x ~ 0], [x])

# specify the harmonic ansatz for x and y: x = u(T) cos(ωt) + v(T) sin(ωt)
add_harmonic!(diff_eq, x, ωₚ)

# implement ansatz to get harmonic equations
harmonic_eq_dummy = get_harmonic_equations(diff_eq)

# Define hilbert space
h = FockSpace(:cavity)

# Define the fundamental operators
@qnumbers a::Destroy(h)

# Define parameters
@variables ωₚ::Real K::Real F::Real κ::Real Δ::Real ω₀::Real # λₑ::Real Δₑ::Real

# Hamiltonian

# H_RWA = (-Δ + 2K) * a' * a + K * (a'^2 * a^2) - G * (a' * a' + a * a) / 2
Gₑ = -(256 * F^2) * K * ω₀^3 / (ωₚ * (Δ^4 - 16 * Δ^2 * ωₚ^2)) # α is taken to be -1
Δₑ = -(512 * F^2) * K * ω₀^3 * (Δ^2 + 16 * ωₚ^2) / (ωₚ * (Δ^3 - 16 * Δ * ωₚ^2)^2) # α is taken to be -1
detuning = (ωₚ^2 - ω₀^2) / (2ωₚ) #- 2 * K * ω₀^2 / ωₚ^2
H_RWA = -(detuning + Δₑ) * a' * a + K * (a'^2 * a^2) - Gₑ * (a' * a' + a * a) / 2

#TODO rescale K
# We do not add the -2K kerr shift as it is included in the bare resonants frequency

# Collapse operators
rates = [κ];
Jop = [a];

# Derive equations
ops = [a, a']

eqs_RWA = meanfield(ops, H_RWA, Jop; rates=rates, order=1)

eqs_completed_RWA = complete(eqs_RWA)

problem_c1 = HarmonicBalance.Problem(eqs_completed_RWA, Num[ωₚ, ω₀, Δ, K, F, κ])
