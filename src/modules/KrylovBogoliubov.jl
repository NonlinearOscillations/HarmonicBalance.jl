module KrylovBogoliubov

using DocStringExtensions
using OrderedCollections: OrderedDict

using Symbolics
using Symbolics:
    unwrap,
    diff2term,
    var_from_nested_derivative,
    lower_varname
using SymbolicUtils: BasicSymbolic, isdiv

using HarmonicBalance
using HarmonicBalance:
    rearrange!,
    flatten,
    is_harmonic,
    _create_harmonic_variable,
    trig_reduce,
    get_independent,
    simplify_complex,
    is_trig,
    substitute_all,
    slow_flow,
    _remove_brackets,
    get_all_terms

include("KrylovBogoliubov/first_order_transform.jl")
include("KrylovBogoliubov/KrylovEquation.jl")

export first_order_transform!,
    is_rearranged_standard, rearrange_standard!, get_equations, get_krylov_equations

end
