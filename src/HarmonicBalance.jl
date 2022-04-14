module HarmonicBalance

    using Printf
    using DataStructures
    import Base: show, display; export show
    export *
    export @variables
    export d
    using Symbolics
    import SymbolicUtils: Term, Add, Div, Mul, Pow, Sym
    using DocStringExtensions

    import Base: ComplexF64, Float64; export ComplexF64, Float64
    ComplexF64(x::Complex{Num}) = ComplexF64(Float64(x.re) + im*Float64(x.im))
    Float64(x::Complex{Num}) = Float64(ComplexF64(x))
    Float64(x::Num) = Float64(x.val)

   # default global settings
   export im_tol
   im_tol = 1E-6
   function set_imaginary_tolerance(x::Float64)
       @eval(im_tol = $x)
   end

   export is_real
    is_real(x) = abs(imag(x)) < im_tol
    is_real(x::Array) = any(is_real.(x))

    # Symbolics does not natively support complex exponentials of variables
    import Base: exp
    exp(x::Complex{Num}) = x.re.val == 0 ? exp(im*x.im.val) : exp(x.re.val + im*x.im.val)

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
    include("plotting_static.jl")
    include("plotting_interactive.jl")
    include("hysteresis_sweep.jl")

    include("modules/HC_wrapper.jl")
    using .HC_wrapper

    include("modules/LinearResponse.jl")
    using .LinearResponse

    include("modules/TimeEvolution.jl")
    using .TimeEvolution

    include("modules/Hopf.jl")
    using .Hopf



end # module
