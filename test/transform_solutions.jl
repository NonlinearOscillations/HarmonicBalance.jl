using HarmonicBalance

@variables Ω γ λ F x θ η α ω0 ω t T ψ
@variables x(t)

natural_equation = d(d(x,t),t) + γ*d(x,t) + Ω^2*x + α*x^3
forces = F*cos(ω*t)
dEOM = DifferentialEquation(natural_equation ~ forces, x)

add_harmonic!(dEOM, x, ω)
harmonic_eq = get_harmonic_equations(dEOM, slow_time=T, fast_time=t);

fixed = (Ω => 1.0, γ => 1E-2, F => 1E-3,  α => 1.)
varied = ω => range(0.9, 1.1, 10)
res = get_steady_states(harmonic_eq, varied, fixed, show_progress=false, seed=SEED);

transform_solutions(res, "u1^2+v1^2")

times = 0:0.1:10
soln = res[5][1]
HarmonicBalance.to_lab_frame(soln, res, times)
HarmonicBalance.to_lab_frame_velocity(soln, res, times)
