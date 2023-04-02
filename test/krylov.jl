using HarmonicBalance

@variables t T x(t) y(t) # symbolic variables
@variables ω ω0 γ F α λ ψ θ η

eq = [d(d(x,t),t) + γ*d(x,t) + ω0^2*(1-λ*cos(2*ω*t+ψ))*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t+θ)]

diff_eom = DifferentialEquation(eq, [x])

add_harmonic!(diff_eom, x, ω) # x will rotate at ω

get_krylov_equations(diff_eom)
