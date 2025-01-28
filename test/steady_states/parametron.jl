using HarmonicBalance
using HarmonicBalance: OrderedDict
using Symbolics
#using Test # do not use Test as this file is used for precompilation

@variables Ω γ λ F x θ η α ω0 ω t T ψ
@variables x(t)

natural_equation =
    d(d(x, t), t) +
    γ * d(x, t) +
    Ω^2 * (1 - λ * cos(2 * ω * t + ψ)) * x +
    α * x^3 +
    η * d(x, t) * x^2
forces = F * cos(ω * t + θ)
dEOM = DifferentialEquation(natural_equation + forces, x)
add_harmonic!(dEOM, x, ω)
harmonic_eq = get_harmonic_equations(dEOM; slow_time=T, fast_time=t);

method = HarmonicBalance.WarmUp(; seed=SEED)

@testset "undriven parametron" begin
    fixed = OrderedDict(
        Ω => 1.0, γ => 1e-2, λ => 5e-2, F => 0, α => 1.0, η => 0.3, θ => 0, ψ => 0
    )
    varied = OrderedDict(ω => range(0.9, 1.1, 20))
    @test substitute(
        sum(harmonic_eq.parameters), merge(Dict(fixed), Dict(varied[ω] => 0))
    ) isa Number

    @testset "Problem" begin
        prob = HarmonicBalance.Problem(harmonic_eq, OrderedDict(varied), OrderedDict(fixed))

        @test length(harmonic_eq.equations) == 2
        @test length(prob.variables) == 2
        @test length(prob.parameters) == 9
        @test length(prob.system.parameters) == 9
        res = get_steady_states(prob, method; show_progress=false)
    end

    @testset "steady states" begin
        res = get_steady_states(harmonic_eq, method, varied, fixed; show_progress=false)
        @test length(res.solutions) == 20
        @test length(res.solutions[1]) == 5
        @test sum(any.(get_class(res, "physical"))) == 5
        @test sum(any.(get_class(res, "stable"))) == 3

        classify_solutions!(res, "sqrt(u1^2 + v1^2) > 1e-6", "nonzero")

        stable_list = get_class(res, "stable")
        zeros_list = get_class(res, "nonzero")

        @test sum(stable_list[1]) == 18
        @test stable_list[2] == stable_list[3]
        @test zeros_list[1] == zeros(Int, 20)
        @test stable_list[2] == stable_list[3]
    end

    @testset "implicit jacobian" begin
        harmonic_eq = get_harmonic_equations(dEOM; jacobian=false)
        p = HarmonicBalance.Problem(harmonic_eq, varied, fixed)
        @test round.(real.(p.jacobian(zeros(3)))) == [-98.0 0.0; 0.0 -102.0]
        res = get_steady_states(p, method; show_progress=false)
    end

    @testset "2d sweep" begin # try to run a 2D calculation
        fixed = (Ω => 1.0, γ => 1e-2, F => 0, α => 1.0, η => 0.1, θ => 0, ψ => 0)
        varied = (ω => range(0.9, 1.1, 20), λ => range(0.01, 0.05, 20))
        res = get_steady_states(harmonic_eq, method, varied, fixed; show_progress=false)
    end
end

@testset "harmonic equation" begin
    @variables u1, v1
    ref1 =
        (Ω^2) * u1 +
        F * cos(θ) +
        γ * Differential(T)(u1) +
        (3//4) * α * (u1^3) +
        γ * ω * v1 +
        (2//1) * ω * Differential(T)(v1) +
        (1//4) * η * ω * (v1^3) +
        (3//4) * η * (u1^2) * Differential(T)(u1) +
        (1//4) * η * (v1^2) * Differential(T)(u1) +
        (3//4) * α * (v1^2) * u1 +
        (1//4) * η * ω * (u1^2) * v1 +
        (1//2) * η * u1 * v1 * Differential(T)(v1) +
        (1//2) * λ * (Ω^2) * v1 * sin(ψ) - (ω^2) * u1 - (1//2) * λ * (Ω^2) * u1 * cos(ψ)
    ref2 =
        γ * Differential(T)(v1) +
        (Ω^2) * v1 +
        (3//4) * α * (v1^3) +
        (3//4) * α * (u1^2) * v1 +
        (1//4) * η * (u1^2) * Differential(T)(v1) +
        (3//4) * η * (v1^2) * Differential(T)(v1) +
        (1//2) * λ * (Ω^2) * v1 * cos(ψ) +
        (1//2) * η * u1 * v1 * Differential(T)(u1) +
        (1//2) * λ * (Ω^2) * u1 * sin(ψ) - F * sin(θ) - (ω^2) * v1 -
        (2//1) * ω * Differential(T)(u1) - γ * ω * u1 - (1//4) * η * ω * (u1^3) -
        (1//4) * η * ω * (v1^2) * u1

    averaged = HarmonicBalance._remove_brackets(harmonic_eq)
    @test isequal(simplify(expand(averaged[1] - ref1)), 0)
    @test isequal(simplify(expand(averaged[2] - ref2)), 0)
end
