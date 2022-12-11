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
    include("hysteresis_sweep.jl")

    include("modules/HC_wrapper.jl")
    using .HC_wrapper

    include("modules/LinearResponse.jl")
    using .LinearResponse

    include("modules/TimeEvolution.jl")
    using .TimeEvolution
    export ParameterSweep, ODEProblem, solve

    include("modules/LimitCycles.jl")
    using .LimitCycles


    @precompile_setup begin
        # Putting some things in `setup` can reduce the size of the
        # precompile file and potentially make loading faster.
        #list = [OtherType("hello"), OtherType("world!")]

        @variables Ω,γ,λ,F, x,θ,η,α, ω0, ω,t,T, ψ
        @variables x(t)

        natural_equation =  d(d(x,t),t) + γ*d(x,t) + Ω^2*(1-λ*cos(2*ω*t+ψ))*x + α * x^3 +η *d(x,t) * x^2
        forces =  F*cos(ω*t+θ)

        @precompile_all_calls begin
            # all calls in this block will be precompiled, regardless of whether
            # they belong to your package or not (on Julia 1.8 and higher)
            dEOM = DifferentialEquation(natural_equation + forces, x)
            add_harmonic!(dEOM, x, ω)
            harmonic_eq = get_harmonic_equations(dEOM);
        end
    end


end # module
