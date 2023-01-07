module HarmonicBalance

    using Printf
    using OrderedCollections
    import Base: show, display; export show
    export *
    export @variables
    export d
    export plot
    using Symbolics
    using ProgressMeter
    import Symbolics.SymbolicUtils: Term, Add, Div, Mul, Pow, Sym
    using DocStringExtensions
    using SnoopPrecompile

    import Base: ComplexF64, Float64; export ComplexF64, Float64
    ComplexF64(x::Complex{Num}) = ComplexF64(Float64(x.re) + im*Float64(x.im))
    Float64(x::Complex{Num}) = Float64(ComplexF64(x))
    Float64(x::Num) = Float64(x.val)

   # default global settings
   export IM_TOL
   IM_TOL::Float64 = 1E-6
   function set_imaginary_tolerance(x::Float64)
       @eval(IM_TOL::Float64 = $x)
   end

   export is_real
    is_real(x) = abs(imag(x)) / abs(real(x)) < IM_TOL::Float64 || abs(x) < 1e-70
    is_real(x::Array) = is_real.(x)

    # Symbolics does not natively support complex exponentials of variables
    import Base: exp
    exp(x::Complex{Num}) = x.re.val == 0 ? exp(im*x.im.val) : exp(x.re.val + im*x.im.val)

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

    include("modules/HC_wrapper.jl")
    using .HC_wrapper

    include("modules/LinearResponse.jl")
    using .LinearResponse

    include("modules/TimeEvolution.jl")
    using .TimeEvolution
    export ParameterSweep, ODEProblem, solve

    include("modules/LimitCycles.jl")
    using .LimitCycles

    precomp_path = (@__DIR__) * "/../test/"
    @precompile_all_calls include(precomp_path * "parametron.jl")
    @precompile_all_calls include(precomp_path * "plotting.jl")

end # module
