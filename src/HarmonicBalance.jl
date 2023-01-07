module HarmonicBalance

    using HarmonicBalanceBase
    # Symbolics_customised
    export quick_cancel
    # Symbolics_utils
    export *
    export @variables
    export rearrange
    export drop_powers
    export get_averaged_equations
    export d
    export substitute_all
    export get_all_terms
    export var_name
    # Types
    export ParameterRange
    export ParameterList
    export StateDict
    export SteadyState
    export ParameterVector
    export DifferentialEquation
    export HarmonicVariable
    export HarmonicEquation
    export Problem
    export Result
    # DE
    export add_harmonic!
    export is_harmonic
    export get_variables
    export get_independent_variables
    export get_harmonic_equations
    export is_rearranged
    export slow_flow, slow_flow!
    export _remove_brackets

    using HarmonicBalanceBase.OrderedCollections
    using HarmonicBalanceBase.Symbolics
    using HarmonicBalanceBase.DocStringExtensions


    import Base: show, display; export show
    import Base: ComplexF64, Float64; export ComplexF64, Float64
    ComplexF64(x::Complex{Num}) = ComplexF64(Float64(x.re) + im*Float64(x.im))
    Float64(x::Complex{Num}) = Float64(ComplexF64(x))
    Float64(x::Num) = Float64(x.val)
    # Symbolics does not natively support complex exponentials of variables
    import Base: exp
    exp(x::Complex{Num}) = x.re.val == 0 ? exp(im*x.im.val) : exp(x.re.val + im*x.im.val)

    using Printf
    using ProgressMeter
    using SnoopPrecompile

   # default global settings
   export IM_TOL
   IM_TOL::Float64 = 1E-6
   function set_imaginary_tolerance(x::Float64)
       @eval(IM_TOL::Float64 = $x)
   end

   export is_real
    is_real(x) = abs(imag(x)) / abs(real(x)) < IM_TOL::Float64 || abs(x) < 1e-70
    is_real(x::Array) = is_real.(x)


    include("types.jl")
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

    # precomp_path = (@__DIR__) * "/../test/"
    # @precompile_all_calls include(precomp_path * "parametron.jl")
    # @precompile_all_calls include(precomp_path * "plotting.jl")

end # module
