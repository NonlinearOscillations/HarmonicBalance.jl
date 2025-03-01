using HarmonicBalance
import HarmonicBalance.LinearResponse.plot_linear_response
using HarmonicBalance: OrderedDict
using Test

@testset "van der Pol oscillator " begin
    @variables ω_lc, t, ω0, x(t), μ

    natural_equation = d(d(x, t), t) - μ * (1 - x^2) * d(x, t) + x
    dEOM = DifferentialEquation(natural_equation, x)

    # order matters for 1*ω_lc gauge to be fixed
    add_harmonic!(dEOM, x, [ω_lc, 3 * ω_lc])

    harmonic_eq = get_harmonic_equations(dEOM)
    HarmonicBalance.LimitCycles._choose_fixed(harmonic_eq, ω_lc)

    fixed = OrderedDict{HarmonicBalance.Num,Float64}()
    varied = OrderedDict(μ => range(2, 3, 2))
    # method = HarmonicBalance.WarmUp(; seed=SEED)
    prob = HarmonicBalance.LimitCycles.limit_cycle_problem(harmonic_eq, varied, fixed, ω_lc)
    result = get_limit_cycles(
        harmonic_eq, WarmUp(), varied, fixed, ω_lc; show_progress=false
    )

    @test sum(any.(get_class(result, "stable"))) == 4
    @test sum(any.(get_class(result, "unique_cycle"))) == 1

    plot(result; y="ω_lc")
    plot_linear_response(result, x, 1; Ω_range=range(2.4, 2.6, 2), order=1)
end
# real(ComplexF64[0.0 + 0.0im 1.375193014698595 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im -0.09353667681429773 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im -0.10927509336716573 + 0.0im 4.158613703389808 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im -4.092544384801762 + 0.0im -0.0777982602614297 + 0.0im])
# takes to long
# @testset "coupled modes" begin
# @variables F, ω, ω_lc, t, x(t), y(t)

# eqs = [d(x,t,2) + 1.0^2*x - x^3 - 0.006*y ~ F*cos(ω*t),
#     d(y,t,2) + 1.005^2*y - y^3 - 0.006*x ~ 0]

# # differential equations
# diffeq = DifferentialEquation(eqs, [x,y])

# # specify the harmonic ansatz for x and y: x = u(T) cos(ωt) + v(T) sin(ωt)
# add_harmonic!(diffeq, x, ω)
# add_harmonic!(diffeq, y, ω)
# add_harmonic!(diffeq, x, ω + ω_lc)
# add_harmonic!(diffeq, y, ω + ω_lc)
# add_harmonic!(diffeq, x, ω - ω_lc)
# add_harmonic!(diffeq, y, ω - ω_lc)

# harmonic_eq = get_harmonic_equations(diffeq);

# fixed = (F => 0.0015)
# varied = (ω => range(0.992, 0.995, 2))

# # results
# result = get_limit_cycles(harmonic_eq, varied, fixed, ω_lc)
# end
