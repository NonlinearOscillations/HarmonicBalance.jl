module LimitCycles

using DocStringExtensions

using Symbolics: Symbolics, Num, expand_derivatives

using HarmonicBalance
using HarmonicBalance:
    HarmonicBalanceMethod,
    WarmUp,
    HomotopyContinuationProblem,
    Result,
    get_steady_states,
    order_branches!,
    find_branch_order,
    _remove_brackets,
    classify_solutions,
    _symidx,
    _is_physical,
    substitute_all,
    var_name,
    get_implicit_Jacobian,
    _free_symbols,
    OrderedDict,
    promote_types,
    JacobianFunction

using HarmonicBalance.ExprUtils: get_all_terms

include("gauge_fixing.jl")
include("analysis.jl")

export get_cycle_variables, get_limit_cycles, add_pairs!

end
