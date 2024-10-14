module LimitCycles

using DocStringExtensions

using Symbolics: Symbolics, Num, expand_derivatives

using HarmonicBalance
using HarmonicBalance:
    order_branches!,
    find_branch_order,
    _remove_brackets,
    classify_solutions,
    _free_symbols,
    _symidx,
    _is_physical,
    substitute_all,
    var_name
using HarmonicBalance.LinearResponse: get_implicit_Jacobian
using HarmonicBalance.ExprUtils: get_all_terms

include("gauge_fixing.jl")
include("analysis.jl")

export get_cycle_variables, get_limit_cycles, add_pairs!

end
