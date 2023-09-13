
using HarmonicBalance
@variables Δω, t, ω0, x(t), μ

natural_equation =  d(d(x,t),t) - μ*(1-x^2) * d(x,t) + x
dEOM = DifferentialEquation(natural_equation, x)

vdp!(k) = add_harmonic!(dEOM, x, (2*k-1)*Δω)
[vdp!(k) for k in 1:2]

harmonic_eq = get_harmonic_equations(dEOM)

fixed = ();
varied = μ => range(1,10,2)

result = HarmonicBalance.LimitCycles.get_steady_states(harmonic_eq, varied, fixed, Δω)

plot(result, y="Δω")

plot_linear_response(result, x, branch=1, Ω_range=range(0.9,1.1,2), order=1)
