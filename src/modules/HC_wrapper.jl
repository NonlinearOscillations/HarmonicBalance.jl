module HC_wrapper

    using HomotopyContinuation: independent_normal
    using Base: get_uuid_name
    using HarmonicBalanceBase.Symbolics
    using HomotopyContinuation
    using HarmonicBalanceBase.DocStringExtensions
    using ..HarmonicBalance

    include("HC_wrapper/homotopy_interface.jl")
end
