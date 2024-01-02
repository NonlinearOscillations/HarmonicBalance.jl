using HarmonicBalance
import HarmonicBalance.LinearResponse.plot_linear_response
import HarmonicBalance.LimitCycles.get_limit_cycles

@variables ω_lc, t, ω0, x(t), μ

natural_equation = d(d(x, t), t) - μ * (1 - x^2) * d(x, t) + x
dEOM = DifferentialEquation(natural_equation, x)

# order should NOT matter if correct gauge is fixed (that corresponding to 1*ω_lc)
add_harmonic!(dEOM, x, [ω_lc, 3 * ω_lc])

harmonic_eq = get_harmonic_equations(dEOM)

fixed = ();
varied = μ => range(1, 5, 5)

result = get_limit_cycles(harmonic_eq, varied, fixed, cycle_harmonic=ω_lc, show_progress=false)

@test sum(any.(classify_branch(result, "stable"))) == 4
@test sum(any.(classify_branch(result, "unique_cycle"))) == 1

plot(result, y="ω_lc")

plot_linear_response(result, x, branch=1, Ω_range=range(0.9, 1.1, 2), order=1)
