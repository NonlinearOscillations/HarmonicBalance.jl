module HarmonicBalance

    using Printf
    using OrderedCollections
    using Symbolics
    using ProgressMeter
    using DocStringExtensions
    using BijectiveHilbert
    using LinearAlgebra
    using Plots, Latexify
    import HomotopyContinuation
    import Distances
    # using SnoopPrecompile

    include("types.jl")

    include("utils.jl")
    include("Symbolics_customised.jl")
    include("Symbolics_utils.jl")
    include("DifferentialEquation.jl")
    include("HarmonicVariable.jl")
    include("HarmonicEquation.jl")
    include("solve_homotopy.jl")
    include("sorting.jl")
    include("classification.jl")
    include("saving.jl")
    include("transform_solutions.jl")
    include("plotting_Plots.jl")

    include("modules/HC_wrapper.jl")
    using .HC_wrapper

    include("modules/LinearResponse.jl")
    using .LinearResponse

    include("modules/LimitCycles.jl")
    using .LimitCycles

    include("modules/KrylovBogoliubov.jl")
    using .KrylovBogoliubov
    export first_order_transform!, is_rearranged_standard, rearrange_standard!, get_equations
    export get_krylov_equations

    using PackageExtensionCompat
    function __init__()
        @require_extensions
    end
    export ParameterSweep, ODEProblem, solve, follow_branch

    # precomp_path = (@__DIR__) * "/../test/"
    # @precompile_all_calls include(precomp_path * "parametron.jl")
    # @precompile_all_calls include(precomp_path * "plotting.jl")

end # module
