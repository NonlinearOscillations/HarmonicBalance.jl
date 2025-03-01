using QuantumCumulants, HarmonicBalance

# first order cumulant

h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
@variables Δ::Real U::Real G::Real κ::Real
param = [Δ, U, G, κ]

H_RWA = (-Δ + U) * a' * a + U * (a'^2 * a^2) / 2 - G * (a' * a' + a * a) / 2
ops = [a, a']

eqs_RWA = meanfield(ops, H_RWA, [a]; rates=[κ], order=1)
eqs_completed_RWA = complete(eqs_RWA)

#

fixed = (U => 0.001, κ => 0.002)
varied = (Δ => range(-0.03, 0.03, 100), G => range(1e-5, 0.02, 100))
problem_c1 = HarmonicBalance.Problem(eqs_completed_RWA, param, varied, fixed)

#

result = get_steady_states(problem_c1, WarmUp())
plot_phase_diagram(result; class="stable")

# second order cumulant

ops = [a]
eqs_RWA = meanfield(ops, H_RWA, Jop; rates=rates, order=2)
eqs_c2 = complete(eqs_RWA)
problem_c2 = HarmonicBalance.Problem(eqs_c2, param, varied, fixed)

#

result = get_steady_states(problem_c2, WarmUp())
plot_phase_diagram(result; class="stable", clim=(0, 4))

#

fixed = (U => 0.001, κ => 0.002, G => 0.01)
varied = (Δ => range(-0.03, 0.03, 200))
problem_c2 = HarmonicBalance.Problem(eqs_c2, param, varied, fixed)
result = get_steady_states(problem_c2, TotalDegree())
plot(result; y="aᵣ")

#

classify_solutions!(result, "a⁺aᵣ < 0", "neg photon number");
plot(result; y="aᵣ", class="stable", not_class="neg photon number")

# third order cumulant

ops = [a]
eqs_RWA = meanfield(ops, H_RWA, Jop; rates=rates, order=3)
eqs_c3 = complete(eqs_RWA)
problem_c3 = HarmonicBalance.Problem(eqs_c3, param, varied, fixed)

#

result = get_steady_states(problem_c3, WarmUp())
plot_phase_diagram(result; class="stable", clim=(0, 4))

#

fixed = (U => 0.001, κ => 0.002, G => 0.01)
varied = (Δ => range(-0.03, 0.03, 50))
problem_c3 = HarmonicBalance.Problem(eqs_c3, param, varied, fixed)
result = get_steady_states(problem_c3, TotalDegree())
plot(result; y="aᵣ")

#

classify_solutions!(result, "a⁺aᵣ < 0", "neg photon number");
plot(result; y="aᵣ", class="stable", not_class="neg photon number")
