using HarmonicBalance, OrderedCollections
using Test, TestExtras

using HarmonicBalance: Problem

@variables α, ω, ω0, F, γ, t, x(t);
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω)
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = OrderedDict(α => 1, ω0 => 1.0, γ => 1e-2, F => 1e-6)
varied = OrderedDict(ω => range(0.9, 1.1, 10))
prob = Problem(harmonic_eq)

@testset "Jacobian FunctionWrapper" begin
    using FunctionWrappers: FunctionWrapper
    using Symbolics: build_function
    using HarmonicBalance: _free_symbols, substitute_all

    complex_vars = [1.0 + 0im, 1.0 + 0im, 1.0 + 0im]
    float_vars = [1.0, 1.0, 1.0]

    J = substitute_all.(prob.jacobian, Ref(fixed))
    jacfunc = build_function(J, _free_symbols(prob, varied); expression=Val(false))[1]
    wrapped_jac = FunctionWrapper{Matrix{ComplexF64},Tuple{Vector{ComplexF64}}}(jacfunc)
    @test wrapped_jac(complex_vars) isa Matrix{ComplexF64}
    @test wrapped_jac(float_vars) isa Matrix{ComplexF64}
    @constinferred wrapped_jac(complex_vars)
    @constinferred wrapped_jac(float_vars)

    jacfunc = build_function(J, _free_symbols(prob, varied); expression=Val(false))[1]
    jacfunc′ = build_function(J, _free_symbols(prob, varied)...; expression=Val(false))[1]
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
