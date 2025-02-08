module HarmonicBalance

# default global settings
IM_TOL::Float64 = 1e-6
function set_imaginary_tolerance(x::Float64)
    @eval(IM_TOL::Float64 = $x)
end

using DocStringExtensions
using JLD2: JLD2
using DelimitedFiles: DelimitedFiles, writedlm
using OrderedCollections: OrderedDict, OrderedSet
using ProgressMeter: ProgressMeter, Progress
using LinearAlgebra: LinearAlgebra, eigvals
using Random: Random # for setting seed
import FunctionWrappers: FunctionWrapper

using Distances: Distances
using BijectiveHilbert: BijectiveHilbert, Simple2D, decode_hilbert!, encode_hilbert
using HomotopyContinuation: HomotopyContinuation
const HC = HomotopyContinuation

using Symbolics:
    Symbolics,
    Num,
    Equation,
    @variables,
    expand_derivatives,
    get_variables,
    Differential,
    unwrap,
    wrap,
    diff2term,
    var_from_nested_derivative,
    lower_varname
using SymbolicUtils: SymbolicUtils

include("modules/ExprUtils/ExprUtils.jl")
using .ExprUtils: is_harmonic, substitute_all, drop_powers, count_derivatives, hasnan

# symbolics equations
export @variables, d
export DifferentialEquation
export HarmonicVariable
export HarmonicEquation

export rearrange_standard
export rearrange_standard!
export first_order_transform!
export is_rearranged_standard
export get_equations

export get_harmonic_equations
export get_krylov_equations
export add_harmonic!

export get_independent_variables
export get_variables

# handle solutions
export get_steady_states
export classify_solutions!
export get_class
export get_single_solution
export transform_solutions
export IM_TOL

# methods
export WarmUp
export TotalDegree
export Polyhedral

# Limit cycles
export get_cycle_variables, get_limit_cycles, add_pairs!

# LinearResponse
export get_Jacobian

# plotting
export plot_linear_response
export plot_phase_diagram
export plot_rotframe_jacobian_response
export plot_eigenvalues
export plot_spaghetti

# extension functions
export AdiabaticSweep
export steady_state_sweep
export plot_1D_solutions_branch
export follow_branch

# src code

include("extension_functions.jl")
include("utils.jl")
include("types.jl")
include("DifferentialEquation.jl")
include("HarmonicVariable.jl")
include("HarmonicEquation.jl")
include("Problem.jl")
include("Jacobian.jl")
include("Result.jl")
include("methods.jl")

include("solve_homotopy.jl")
include("sorting.jl")
include("classification.jl")

include("saving.jl")
include("transform_solutions.jl")
# include("plotting.jl")

include("modules/HC_wrapper.jl")
using .HC_wrapper

include("modules/LinearResponse/LinearResponse.jl")
using .LinearResponse

include("modules/LimitCycles/LimitCycles.jl")
using .LimitCycles

include("modules/KrylovBogoliubov.jl")
using .KrylovBogoliubov

# Precompilation setup
using PrecompileTools: @setup_workload, @compile_workload
# @setup_workload begin
#     # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
#     # precompile file and potentially make loading faster.
#     @compile_workload begin
#         # all calls in this block will be precompiled, regardless of whether
#         # they belong to your package or not (on Julia 1.8 and higher)
#         include("precompilation.jl")
#     end
# end

# Error hint for extensions stubs
function __init__()
    Base.Experimental.register_error_hint(
        _error_hinter("OrdinaryDiffEq", :TimeEvolution, follow_branch), MethodError
    )
    Base.Experimental.register_error_hint(
        _error_hinter("OrdinaryDiffEq", :TimeEvolution, plot_1D_solutions_branch),
        MethodError,
    )
    Base.Experimental.register_error_hint(
        _error_hinter("SteadyStateDiffEq", :SteadyStateDiffEqExt, steady_state_sweep),
        MethodError,
    )
    for func in [
        plot_spaghetti,
        plot_eigenvalues,
        plot_rotframe_jacobian_response,
        plot_phase_diagram,
        plot_linear_response,
    ]
        Base.Experimental.register_error_hint(
            _error_hinter("Plots", :PlotsExt, func), MethodError
        )
    end
    return nothing
end

end # module
