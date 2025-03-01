using HarmonicBalance, BenchmarkTools
using Test

using HarmonicBalance: Problem, substitute_all, _free_symbols, OrderedDict
using Symbolics: build_function
using FunctionWrappers: FunctionWrapper

@variables α, ω, ω0, F, γ, t, x(t);
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω)
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = OrderedDict(α => 1, ω0 => 1.0, γ => 1e-2, F => 1e-6)
varied = OrderedDict(ω => range(0.9, 1.1, 10))
prob = Problem(harmonic_eq)

J = substitute_all.(prob.jacobian, Ref(fixed))
jacfunc = build_function(J, _free_symbols(prob); expression=Val(false))[1]
jacfunc′ = build_function(J, _free_symbols(prob)...; expression=Val(false))[1]

struct JacobianFunctionTest{T}
    J::FunctionWrapper{Matrix{T},Tuple{Vector{T}}}
end
(cb::JacobianFunctionTest)(v) = cb.J(v)

struct JacobianFunctionTest′{T,N}
    J::FunctionWrapper{Matrix{T},NTuple{N,T}}
end
(cb::JacobianFunctionTest′)(v...) = cb.J(v...)

wrapped_jac = JacobianFunctionTest{ComplexF64}(jacfunc)
wrapped_jac′ = JacobianFunctionTest′{ComplexF64,3}(jacfunc′)

data = rand(ComplexF64, 3)
@btime wrapped_jac($data) # 170.751 ns (5 allocations: 256 bytes)
@btime jacfunc($data) # 163.393 ns (5 allocations: 256 bytes)

a = rand(ComplexF64)
b = rand(ComplexF64)
c = rand(ComplexF64)
@btime wrapped_jac′($a, $b, $c) # 310.976 ns (8 allocations: 352 bytes)
@btime jacfunc′($a, $b, $c) # 159.296 ns (8 allocations: 352 bytes)

@variables α, ω, ω0, F, γ, t, J1, x(t), y(t), z(t);
diff_eq = DifferentialEquation(
    [
        d(x, t, 2) + ω0 * x + α * x^3 + γ * d(x, t) + J1 * y + J1 * z ~ F * cos(ω * t),
        d(y, t, 2) + ω0 * y + α * y^3 + γ * d(y, t) + J1 * x + J1 * z ~ F * cos(ω * t),
        d(z, t, 2) + ω0 * z + α * z^3 + γ * d(z, t) + J1 * x + J1 * y ~ F * cos(ω * t),
    ],
    [x, y, z],
)
add_harmonic!(diff_eq, x, ω)
add_harmonic!(diff_eq, y, ω)
add_harmonic!(diff_eq, z, ω)

harmonic_eq = get_harmonic_equations(diff_eq)
prob = Problem(harmonic_eq)

fixed = OrderedDict(α => 1, ω0 => 1.0, γ => 1e-2, F => 1e-6, J1 => 1e-2)
J = substitute_all.(prob.jacobian, Ref(fixed))
jacfunc = build_function(J, _free_symbols(prob); expression=Val(false))[1]
jacfunc′ = build_function(J, _free_symbols(prob)...; expression=Val(false))[1]

wrapped_jac = JacobianFunctionTest{ComplexF64}(jacfunc)
wrapped_jac′ = JacobianFunctionTest′{ComplexF64,7}(jacfunc′)

data = rand(ComplexF64, 7)
@btime wrapped_jac($data) # 14.900 μs (54 allocations: 3.03 KiB)
@btime jacfunc($data) # 14.900 μs (54 allocations: 3.03 KiB)

a = rand(ComplexF64)
b = rand(ComplexF64)
c = rand(ComplexF64)
d1 = rand(ComplexF64)
e = rand(ComplexF64)
f = rand(ComplexF64)
g = rand(ComplexF64)
@btime wrapped_jac′($a, $b, $c, $d1, $e, $f, $g) # 15.100 μs (61 allocations: 3.25 KiB)
@btime jacfunc′($a, $b, $c, $d1, $e, $f, $g) # 14.800 μs (61 allocations: 3.25 KiB)
