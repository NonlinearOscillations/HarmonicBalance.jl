using HarmonicBalance, Symbolics; HB = HarmonicBalance
using Symbolics.Rewriters
using Symbolics: BasicSymbolic, unwrap, arguments, isdiv, isadd, ismul
import Symbolics.SymbolicUtils: @compactified

@variables t T x(t) y(t) # symbolic variables
@variables ω ω0 γ F α λ ψ θ η

eq = [d(d(x,t),t) + γ*d(x,t) + ω0^2*(1-λ*cos(2*ω*t))*x + α*x^3 #= + η*d(x,t)*x^2 ~ F*cos(ω*t+θ)=#]

diff_eom = DifferentialEquation(eq, [x])

add_harmonic!(diff_eom, x, ω) # x will rotate at ω

harmonic_eq = get_krylov_equations(diff_eom, order = 2)

fixed = (ω0 => 1.0, γ => 0.005, α => 1.0, η => 0, F => 0.0, ψ => 0.0, θ => 0.0, λ => 0.02)
varied = (ω => range(0.99, 1.01, 50))
res = get_steady_states(harmonic_eq, varied, fixed, threading =true)

# plot 1D result
plot(res, x = "ω", y = "√(u1^2+v1^2)", class="stable")

fixed = (ω0 => 1.0, γ => 0.005, α => 1.0, η => 0, F => 0.0, ψ => 0.0, θ => 0.0)
varied = (ω => range(0.99, 1.01, 50), λ => range(1e-6, 0.05, 50))
res = get_steady_states(harmonic_eq, varied, fixed, threading =true)

# plot 2D result
plot_phase_diagram(res, class="stable")
