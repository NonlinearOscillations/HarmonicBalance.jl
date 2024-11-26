using QuantumCumulants, HarmonicBalance
using SteadyStateDiffEq, ModelingToolkit
using Plots, LinearAlgebra

# Define hilbert space
h = FockSpace(:cavity)

# Define the fundamental operators
@qnumbers a::Destroy(h)

# Define parameters
@variables Δ::Real U::Real F::Real G::Real κ::Real

# Hamiltonian

H_RWA = (-Δ + U) * a' * a + U * (a'^2 * a^2) / 2 - G * (a' * a' + a * a) / 2

# Collapse operators
rates = [κ];
Jop = [a];

# Derive equations
ops = [a, a']

eqs_RWA = meanfield(ops, H_RWA, Jop; rates=rates, order=1)

eqs_completed_RWA = complete(eqs_RWA)

problem_c1 = HarmonicBalance.Problem(eqs_completed_RWA, Num[Δ, U, G, κ])
hasproperty(problem_c1, :eom)


problem_c1.eom
problem_c1.eom

fixed = (U => 0.001, κ => 0.002, G => 0.01)
varied = (Δ => range(-0.03, 0.03, 50))

# results
result = get_steady_states(problem_c1, WarmUp(), varied, fixed)
