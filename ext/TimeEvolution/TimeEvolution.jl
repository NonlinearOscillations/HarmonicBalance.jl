module TimeEvolution

using DocStringExtensions
using Symbolics: Num, substitute, unwrap
using OrdinaryDiffEqTsit5: OrdinaryDiffEqTsit5

using HarmonicBalance:
    HarmonicBalance,
    StateDict,
    HarmonicEquation,
    _apply_mask,
    _get_mask,
    rearrange_standard,
    is_rearranged,
    filter_duplicate_parameters,
    _parse_expression,
    _set_Plots_default,
    Result,
    substitute_all,
    get_variables,
    transform_solutions,
    get_single_solution,
    follow_branch
const HB = HarmonicBalance

include("sweeps.jl")
include("ODEProblem.jl")
include("hysteresis_sweep.jl")
include("plotting.jl")

export ParameterSweep
export transform_solutions, plot, plot!
export ODEProblem, solve
export follow_branch

end
