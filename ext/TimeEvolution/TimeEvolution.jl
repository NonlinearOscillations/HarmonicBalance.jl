module TimeEvolution

using DocStringExtensions
using Symbolics: Num, substitute, unwrap
using OrdinaryDiffEq: OrdinaryDiffEq

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

export FFT
export ParameterSweep
export transform_solutions, plot, plot!, is_stable
export ODEProblem, solve
export plot_1D_solutions_branch, follow_branch

end
