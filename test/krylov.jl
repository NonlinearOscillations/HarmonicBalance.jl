using HarmonicBalance

@variables t T x(t) y(t) # symbolic variables
@variables ω ω0 γ F

eq = [d(d(x, t),t) + ω0^2*x ~ - γ*d(x,t)]

diff_eom = DifferentialEquation(eq, [x, y])

add_harmonic!(diff_eom, x, ω) # x will rotate at ω


get_krylov_equations(diff_eom)
