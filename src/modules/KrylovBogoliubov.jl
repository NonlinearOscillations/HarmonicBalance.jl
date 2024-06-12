module KrylovBogoliubov

using ..HarmonicBalance
using Symbolics
using LinearAlgebra
using OrderedCollections
using DocStringExtensions

using Symbolics:
    unwrap,
    operation,
    arguments,
    issym,
    diff2term,
    isdiv,
    BasicSymbolic,
    var_from_nested_derivative,
    lower_varname
using HarmonicBalance: is_rearranged, rearrange!, rearrange
using HarmonicBalance: flatten, is_harmonic, _create_harmonic_variable, HarmonicEquation
using HarmonicBalance: trig_reduce, get_independent, simplify_complex, get_Jacobian, is_trig

include("KrylovBogoliubov/first_order_transform.jl")
include("KrylovBogoliubov/KrylovEquation.jl")
end
