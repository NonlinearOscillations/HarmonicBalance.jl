using HarmonicBalance, Symbolics
HB = HarmonicBalance

@variables α, ω, ω0, F, γ, t, x(t);

diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω)
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = (α => 1, ω0 => 1.0, γ => 1e-2, F => 1e-6)
varied = ω => range(0.9, 1.1, 10)

result = get_steady_states(harmonic_eq, varied, fixed; show_progress=false)

@testset "first order" begin
    plot_linear_response(
        result, x; branch=1, Ω_range=range(0.9, 1.1, 10), order=1, logscale=true
    )
    plot_rotframe_jacobian_response(
        result; Ω_range=range(0.01, 1.1, 10), branch=1, logscale=true
    )
end

@testset "second order" begin
    response_matrix = HB.LinearResponse.ResponseMatrix(result)
    M = response_matrix.matrix
    @test M[1](ones(4)) isa ComplexF64

    plot_linear_response(
        result, x; branch=1, Ω_range=range(0.9, 1.1, 10), order=2, logscale=true
    )
end

@testset "eigenvalues" begin
    plot_eigenvalues(result; branch=1)
    plot_eigenvalues(result; branch=1, type=:re, class="all")

    @testset "NaN Exception Error" begin
        @variables α λ ω0 ω ωₚ F t x(t)
        diff_eq = DifferentialEquation(
            d(x, t, 2) + (ω0^2 - λ * cos(2 * ω * t)) * x + α * x^3 + γ * d(x, t) ~
                F * cos(ωₚ * t),
            x,
        )

        add_harmonic!(diff_eq, x, ω)
        add_harmonic!(diff_eq, x, ωₚ)
        harmonic_eq = get_harmonic_equations(diff_eq)

        fixed = (γ => 0.008, ω0 => 1.0, α => 1.0, F => 0.0, ωₚ => 1.0, λ => 0.016)
        varied = (ω => range(0.995, 1.005, 20))
        result_ω = get_steady_states(harmonic_eq, varied, fixed; show_progress=false)
        @test_throws ErrorException plot_eigenvalues(result_ω; branch=1)
    end
end
