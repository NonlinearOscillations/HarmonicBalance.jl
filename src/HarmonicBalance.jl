module HarmonicBalance

using DocStringExtensions
using JLD2: JLD2
using DelimitedFiles: DelimitedFiles, writedlm
using OrderedCollections: OrderedDict, OrderedSet
using ProgressMeter: ProgressMeter, Progress
using LinearAlgebra: eigvals
using Random: Random # for setting seed

using Distances: Distances
using BijectiveHilbert: BijectiveHilbert, Simple2D, decode_hilbert!, encode_hilbert
using HomotopyContinuation: HomotopyContinuation
const HC = HomotopyContinuation

using Plots: Plots, plot, plot!, savefig, heatmap, Plot
using Latexify: Latexify, latexify

using PrecompileTools: @setup_workload, @compile_workload

# default global settings
IM_TOL::Float64 = 1E-6
function set_imaginary_tolerance(x::Float64)
    @eval(IM_TOL::Float64 = $x)
end

include("Symbolics_customised.jl")
include("Symbolics_utils.jl")

include("modules/extention_functions.jl")
include("utils.jl")
include("types.jl")

include("DifferentialEquation.jl")
include("HarmonicVariable.jl")
include("HarmonicEquation.jl")
include("solve_homotopy.jl")
include("sorting.jl")
include("classification.jl")
include("saving.jl")
include("transform_solutions.jl")
include("plotting_Plots.jl")

export show, *, @variables, d, ComplexF64, Float64, IM_TOL

export ParameterRange, ParameterList, StateDict, SteadyState, ParameterVector
export DifferentialEquation, HarmonicVariable, HarmonicEquation, Problem, Result
export get_steady_states, get_single_solution, get_harmonic_equations, add_harmonic!
export get_variables, get_independent_variables, classify_branch, classify_solutions!
export rearrange_standard

export plot, plot!, plot_phase_diagram, savefig, plot_spaghetti

export ParameterSweep, steady_state_sweep
export plot_1D_solutions_branch, follow_branch

include("modules/HC_wrapper.jl")
using .HC_wrapper

include("modules/LinearResponse.jl")
using .LinearResponse
export plot_linear_response, plot_rotframe_jacobian_response, get_Jacobian, plot_eigenvalues
export transform_solutions

include("modules/LimitCycles.jl")
using .LimitCycles
export get_cycle_variables, get_limit_cycles, add_pairs!

include("modules/KrylovBogoliubov.jl")
using .KrylovBogoliubov
export first_order_transform!, is_rearranged_standard, rearrange_standard!, get_equations
export get_krylov_equations

include("modules/FFTWExt.jl")
using .FFTWExt

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        include("precompilation.jl")
    end
end

end # module
