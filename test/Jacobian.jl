using HarmonicBalance, OrderedCollections
using Test, TestExtras

using HarmonicBalance: Problem

@variables α, ω, ω0, F, γ, t, x(t);
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω)
eom = get_harmonic_equations(diff_eq)

fixed = OrderedDict(α => 1, ω0 => 1.0, γ => 1e-2, F => 1e-6)
varied = OrderedDict(ω => range(0.9, 1.1, 10))
prob = Problem(eom, varied, fixed)

@testset "Jacobian FunctionWrapper" begin
    using FunctionWrappers: FunctionWrapper
    using Symbolics: build_function
    using HarmonicBalance: _free_symbols, substitute_all

    complex_vars = [1.0 + 0im, 1.0 + 0im, 1.0 + 0im]
    float_vars = [1.0, 1.0, 1.0]

    J = substitute_all.(eom.jacobian, Ref(fixed))
    jacfunc = build_function(J, _free_symbols(prob); expression=Val(false))[1]
    wrapped_jac = FunctionWrapper{Matrix{ComplexF64},Tuple{Vector{ComplexF64}}}(jacfunc)
    @test wrapped_jac(complex_vars) isa Matrix{ComplexF64}
    @test wrapped_jac(float_vars) isa Matrix{ComplexF64}
    @constinferred wrapped_jac(complex_vars)
    @constinferred wrapped_jac(float_vars)

    jacfunc = build_function(J, _free_symbols(prob); expression=Val(false))[1]
    jacfunc′ = build_function(J, _free_symbols(prob)...; expression=Val(false))[1]
    @test jacfunc(complex_vars) isa Matrix{ComplexF64}
    @test jacfunc(float_vars) isa Matrix{Float64}
    @constinferred jacfunc(complex_vars)
    @constinferred jacfunc(float_vars)
    @test jacfunc′(complex_vars...) isa Matrix{ComplexF64}
    @test jacfunc′(float_vars...) isa Matrix{Float64}
    @constinferred jacfunc′(complex_vars...)
    @constinferred jacfunc′(float_vars...)

    struct JacobianFunctionTest{T}
        J::FunctionWrapper{Matrix{T},Tuple{Vector{T}}}
    end
    (cb::JacobianFunctionTest)(v) = cb.J(v)

    struct JacobianFunctionTest′{T,N}
        J::FunctionWrapper{Matrix{T},NTuple{N,T}}
    end
    (cb::JacobianFunctionTest′)(v...) = cb.J(v...)

    wrapped_jac = JacobianFunctionTest{ComplexF64}(jacfunc)
    # wrapped_jac′ = JacobianFunctionTest′{ComplexF64,3}((x...) -> jacfunc′(x...))
    wrapped_jac′ = JacobianFunctionTest′{ComplexF64,3}(jacfunc′)
    @test wrapped_jac(complex_vars) isa Matrix{ComplexF64}
    @test wrapped_jac(float_vars) isa Matrix{ComplexF64}
    @constinferred wrapped_jac(complex_vars)
    @constinferred wrapped_jac(float_vars)
    @test wrapped_jac′(complex_vars...) isa Matrix{ComplexF64}
    @test wrapped_jac′(float_vars...) isa Matrix{ComplexF64}
end

@testset "NaNMath" begin
    # https://github.com/QuantumEngineeredSystems/HarmonicBalance.jl/issues/357
    using HarmonicBalance

    @variables t x1(t) x2(t)
    @variables b c m γ ψ

    equations = [
        m * d(d(x1, t), t) +
        2 * b * x1 +
        γ * d(x1, t) +
        4 * x1^3 +
        3 * c * Base.cos(ψ)^3 * x1^2 +
        4 * x2^2 * x1 +
        3 * c * (Base.cos(ψ) * Base.sin(ψ)^2 - 1) * x2^2 -
        6 * c * Base.cos(ψ)^2 * Base.sin(ψ) * x1 * x2,
        m * d(d(x2, t), t) + 2 * b * x2 + γ * d(x2, t) + 4 * x2^3 -
        3 * c * Base.sin(ψ)^3 * x2^2 +
        c * (-6 + 3 / 2 * Base.cos(ψ) - 3 / 2 * Base.cos(3 * ψ)) * x1 * x2 +
        4 * x1^2 * x2 - 3 * c * Base.cos(ψ)^2 * Base.sin(ψ) * x1^2,
    ]
    system = DifferentialEquation(equations, [x1, x2])

    add_harmonic!(system, x1, 0)
    add_harmonic!(system, x2, 0)

    harmonic_normal = get_harmonic_equations(system)

    fixed = (γ => 0.01, m => 1, c => -3.5, b => -2) # fixed parameters
    varied = (ψ => range(0.01, 2, 100)) # range of parameter values

    method = TotalDegree()
    result_asym = get_steady_states(
        harmonic_normal, method, varied, fixed; show_progress=false
    )

    @testset "smaller test" begin
        using Symbolics
        m = [cos(ψ) -sin(ψ); sin(ψ) cos(ψ)]
        jacfunc = Symbolics.build_function(m, ψ; expression=Val(false))[1]
        jacfunc((0 + 0.1im))
    end
end
