using QuantumCumulants, HarmonicBalance
using Plots

@testset "KPO" begin
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)
    @variables Δ::Real U::Real F::Real G::Real κ::Real

    H_RWA = (-Δ + U) * a' * a + U * (a'^2 * a^2) / 2 - G * (a' * a' + a * a) / 2
    ops = [a, a']

    eqs_RWA = meanfield(ops, H_RWA, [a]; rates=[κ], order=1)
    eqs_completed_RWA = complete(eqs_RWA)

    fixed = (U => 0.001, κ => 0.00, Δ => 0.0)
    varied = (G => range(0.01, 0.02, 10))
    problem = HarmonicBalance.Problem(eqs_completed_RWA, [Δ, U, G, κ], varied, fixed)
    result = get_steady_states(problem, TotalDegree())
    @test sum(all.(get_class(result, "stable"))) == 2

    fixed = (U => 0.001, κ => 0.00, G => 0.01)
    varied = (Δ => range(-0.03, 0.03, 10))
    problem = HarmonicBalance.Problem(eqs_completed_RWA, [Δ, U, G, κ], varied, fixed)
    result = get_steady_states(problem, TotalDegree())
    @test sum(any.(get_class(result, "stable"))) == 3
end
