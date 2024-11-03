module KrylovBogoliubov

using DocStringExtensions
using OrderedCollections: OrderedDict

using Symbolics
using Symbolics: unwrap, diff2term, var_from_nested_derivative, lower_varname
using SymbolicUtils: BasicSymbolic, isdiv

using HarmonicBalance
using HarmonicBalance:
    rearrange!,
    flatten,
    _create_harmonic_variable,
    slow_flow,
    _remove_brackets,
    get_variables_nums

using HarmonicBalance.ExprUtils:
    get_all_terms,
    substitute_all,
    trig_reduce,
    get_independent,
    simplify_complex,
    is_trig,
    is_harmonic
include("first_order_transform.jl")
include("KrylovEquation.jl")

export first_order_transform!,
    is_rearranged_standard, rearrange_standard!, get_equations, get_krylov_equations

end
