using HarmonicBalance, Symbolics
using Test
HB = HarmonicBalance

@testset "Lorentzian" begin
    using HarmonicBalance.LinearResponse: Lorentzian
    Γ = Lorentzian(1.0, 0.005, 1.0)
    @test typeof(Γ).parameters[1] == Float64
    Γr′ = 2.0 * Γ
    Γl′ = Γ * 2.0
    @test Γr′.A == 2.0 && Γl′.A == 2.0
end

@testset "JacobianSpectrum" begin
    using HarmonicBalance.LinearResponse: JacobianSpectrum, add_peak, evaluate
    s = JacobianSpectrum{Float64}()
    p = Lorentzian(1.0, 0.005, 1.0)
    s = add_peak(s, p)
    @test length(s.peaks) == 1
    s = add_peak(s, p)
    @test length(s.peaks) == 2
    @test (2.0 * s).peaks[1].A == 2.0
    s2 = JacobianSpectrum{Float64}()
    s2 = add_peak(s2, p)
    s = add_peak(s, s2)
    @test length(s.peaks) == 3
    @test evaluate(s, 1.0) == 600.0
end

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
    plot_linear_response(result, x, 1; Ω_range=range(0.9, 1.1, 10), order=1, logscale=true)
    plot_rotframe_jacobian_response(result, 1; Ω_range=range(0.01, 1.1, 10), logscale=true)
end

@testset "second order" begin
    response_matrix = HB.LinearResponse.ResponseMatrix(result)
    M = response_matrix.matrix
    @test M[1](ones(4)) isa ComplexF64

    plot_linear_response(result, x, 1; Ω_range=range(0.9, 1.1, 10), order=2, logscale=true)
end

@testset "second order krylov" begin
    @variables α, ω, ω0, F, γ, t, x(t)
    diff_eq = DifferentialEquation(d(x, t, 2) + ω0^2 * x + α * x^3 ~ F * cos(ω * t), x)
    add_harmonic!(diff_eq, x, ω)
    kylov_eq = get_krylov_equations(diff_eq; order=1)

    fixed = (α => 1.0, ω0 => 1.0, F => 0.002)
    varied = ω => range(0.95, 1.1, 10)
    result = get_steady_states(kylov_eq, varied, fixed)

    Ω_range = range(0.95, 1.1, 10)
    HarmonicBalance.get_linear_response(result, x, Ω_range, 1; order=2)
end

@testset "eigenvalues" begin
    plot_eigenvalues(result, 1)
    plot_eigenvalues(result, 1; type=:re, class="all")

    @testset "NaN Exception Error" begin
        @variables α λ ω0 ω ωₚ γ F t x(t)
        diff_eq = DifferentialEquation(
            d(x, t, 2) + (ω0^2 - λ * cos(2 * ω * t)) * x + α * x^3 + γ * d(x, t), x
        )

        add_harmonic!(diff_eq, x, ω)
        harmonic_eq = get_harmonic_equations(diff_eq)

        fixed = (γ => 0.002, ω0 => 1.0, α => 1.0, λ => 0.016)
        varied = (ω => range(0.995, 1.005, 20))
        result_ω = get_steady_states(harmonic_eq, varied, fixed)
        plot(result_ω; y="u1")
        broken = true
        @test try
            plot_eigenvalues(result_ω; branch=1)
            false
        catch e
            typeof(e) == ErrorException
        end broken = true

        # replace with system where branch has non-physical solution. I think coupled parametrons has this
    end
end
