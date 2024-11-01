module HarmonicBalance

# default global settings
IM_TOL::Float64 = 1E-6
function set_imaginary_tolerance(x::Float64)
    @eval(IM_TOL::Float64 = $x)
end

using DocStringExtensions
using JLD2: JLD2
using DelimitedFiles: DelimitedFiles, writedlm
using OrderedCollections: OrderedDict, OrderedSet
using ProgressMeter: ProgressMeter, Progress
using LinearAlgebra: eigvals
using Random: Random # for setting seed
using EndpointRanges: EndpointRanges

using Distances: Distances
using BijectiveHilbert: BijectiveHilbert, Simple2D, decode_hilbert!, encode_hilbert
using HomotopyContinuation: HomotopyContinuation
const HC = HomotopyContinuation

using Plots: Plots, plot, plot!, savefig, heatmap, Plot
using Latexify: Latexify, latexify

using Symbolics:
    Symbolics,
    Num,
    Equation,
    @variables,
    expand_derivatives,
    get_variables,
    Differential,
    unwrap,
    wrap
using SymbolicUtils: SymbolicUtils

include("ExprUtils/ExprUtils.jl")
using .ExprUtils: is_harmonic, substitute_all, drop_powers, count_derivatives

include("extention_functions.jl")
include("utils.jl")
include("types.jl")
include("DifferentialEquation.jl")
include("HarmonicVariable.jl")
include("HarmonicEquation.jl")
include("Problem.jl")
include("Result.jl")
include("methods.jl")

include("solve_homotopy.jl")
include("sorting.jl")
include("classification.jl")

include("saving.jl")
include("transform_solutions.jl")
include("plotting_Plots.jl")

export show, *, @variables, d, ComplexF64, Float64, IM_TOL

export DifferentialEquation, HarmonicVariable, HarmonicEquation
export get_steady_states, get_single_solution, get_harmonic_equations, add_harmonic!
export get_variables, get_independent_variables, classify_branch, classify_solutions!
export rearrange_standard

export plot, plot!, plot_phase_diagram, savefig, plot_spaghetti

export AdiabaticSweep, steady_state_sweep
export plot_1D_solutions_branch, follow_branch

include("HC_wrapper/HC_wrapper.jl")
using .HC_wrapper

include("LinearResponse/LinearResponse.jl")
using .LinearResponse
export plot_linear_response, plot_rotframe_jacobian_response, get_Jacobian, plot_eigenvalues
export transform_solutions

include("LimitCycles/LimitCycles.jl")
using .LimitCycles
export get_cycle_variables, get_limit_cycles, add_pairs!

include("KrylovBogoliubov/KrylovBogoliubov.jl")
using .KrylovBogoliubov
export first_order_transform!, is_rearranged_standard, rearrange_standard!, get_equations
export get_krylov_equations

include("FFTWExt.jl")
using .FFTWExt

using PrecompileTools: @setup_workload, @compile_workload

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
