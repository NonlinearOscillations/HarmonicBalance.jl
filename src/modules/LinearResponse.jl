module LinearResponse

    using Symbolics: variables
    using LinearAlgebra
    using Printf
    using Symbolics
    using OrderedCollections
    using ..HarmonicBalance
    using ..HC_wrapper
    using DocStringExtensions

    import Base: show; export show

    include("LinearResponse/types.jl")
    include("LinearResponse/utils.jl")
    include("LinearResponse/jacobians.jl")
    include("LinearResponse/Lorentzian_spectrum.jl")
    include("LinearResponse/response.jl")
    include("LinearResponse/plotting.jl")

export get_Jacobian
export plot_linear_response, plot_rotframe_jacobian_response, plot_eigenvalues

end
