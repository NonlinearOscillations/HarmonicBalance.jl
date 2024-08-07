using BenchmarkTools
using HarmonicBalance;
HB = HarmonicBalance;
# using Plots

@variables α Δ ω0 δ F γ η t x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) + η * x^2 * d(x, t) ~
        F * (cos((ω0 + δ - Δ / 2) * t) + cos((ω0 + δ + Δ / 2) * t)),
    x,
)

add_harmonic!(diff_eq, x, ω0 + δ - Δ / 2)
add_harmonic!(diff_eq, x, ω0 + δ + Δ/2)
# add_harmonic!(diff_eq, x, ω0 + δ)

resolution = 2
δ_range = range(-0.05 * 1.25, 0.01 * 1.25, 10);
fixed = (Δ => 1.25, γ => 0.002, ω0 => 41.140, α => -3.2, η => 0.0078, F => 20)
varied_δF = (δ => δ_range)

harmonic_eq = get_harmonic_equations(diff_eq)
prob = HB.Problem(harmonic_eq)

swept_parameters = HB.ParameterRange(varied_δF)
fixed_parameters = HB.ParameterList(fixed)
unique_fixed = HB.filter_duplicate_parameters(swept_parameters, fixed_parameters)
J = HB._compile_Jacobian(prob, swept_parameters, unique_fixed)
# @btime HB._compile_Jacobian($prob, $swept_parameters, $unique_fixed) # 71 ms
result_δF = get_steady_states(harmonic_eq, varied_δF, fixed);
@btime result_δF = get_steady_states($harmonic_eq, $varied_δF, $fixed);
