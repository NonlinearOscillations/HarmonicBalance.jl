using HarmonicBalance
import HarmonicBalance.LinearResponse.plot_linear_response
import HarmonicBalance.LimitCycles.get_limit_cycles


@testset "van der Pol oscillator " begin
    @variables ω_lc, t, ω0, x(t), μ

    natural_equation = d(d(x, t), t) - μ * (1 - x^2) * d(x, t) + x
    dEOM = DifferentialEquation(natural_equation, x)

    # order matters for 1*ω_lc gauge to be fixed
    add_harmonic!(dEOM, x, [ω_lc, 3 * ω_lc])

    harmonic_eq = get_harmonic_equations(dEOM)
    HarmonicBalance.LimitCycles._choose_fixed(harmonic_eq, ω_lc)

    fixed = ();
    varied = μ => range(1, 5, 5)

    result = get_limit_cycles(harmonic_eq, varied, fixed, ω_lc; show_progress=false, seed=SEED)

    @test sum(any.(classify_branch(result, "stable"))) == 4
    @test sum(any.(classify_branch(result, "unique_cycle"))) == 1

    plot(result, y="ω_lc")
    plot_linear_response(result, x, branch=1, Ω_range=range(0.9, 1.1, 2), order=1)
end
