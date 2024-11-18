using HarmonicBalance, OrderedCollections, Symbolics
using Test

@variables α, ω, ω0, F, γ, t, x(t);

diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω)
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = OrderedDict(α => 1, ω0 => 1.0, γ => 1e-2, F => 1e-6)
varied = OrderedDict(ω => range(0.9, 1.1, 10))
prob = HarmonicBalance.Problem(harmonic_eq)

solutions = Vector{Vector{Vector{Complex{Float64}}}}(undef, 10)
seed = 0xd8e5d8df

J = HarmonicBalance.substitute_all.(prob.jacobian, Ref(fixed))
jacfunc = Symbolics.build_function(
    J, HarmonicBalance._free_symbols(prob, varied); expression=Val(false)
)[1]
jacfunc(ComplexF64.([1.0, 1.0, 1.0]))
wrapped_jac =  HarmonicBalance.JacobianFunction{Matrix{ComplexF64}, Tuple{Vector{ComplexF64}}}(x -> jacfunc(x))

@test rapped_jac([1.0, 1.0, 1.0].+0im) isa Matrix{ComplexF64}
@test wrapped_jac([1.0, 1.0, 1.0]) isa Matrix{ComplexF64}

# HarmonicBalance.Result(
#     solutions, varied, fixed, prob, Dict(), zeros(Int64, size(solutions)...), wrapped_jac, seed
# )
# g = FunctionWrapper{Matrix{Float64}, NTuple{3,Float64}}(x -> jacfunc(x...))

# using FunctionWrappers
# import FunctionWrappers: FunctionWrapper

# # For a function that sends (x1::T1, x2::T2, ...) -> ::TN, you use
# # a FunctionWrapper{TN, Tuple{T1, T2, ...}}.
# struct TypeStableStruct
#   fun::FunctionWrapper{Float64, Tuple{Float64, Float64}}
#   second_arg::Float64
# end

# evaluate_strfun(str, arg) = str.fun(arg, str.second_arg)

# example = TypeStableStruct(hypot, 1.0)

# @code_warntype evaluate_strfun(example, 1.5) # all good
