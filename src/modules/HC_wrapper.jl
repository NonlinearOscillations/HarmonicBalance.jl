module HC_wrapper

    using HomotopyContinuation: independent_normal
    using Base: get_uuid_name
    using Symbolics
    using HomotopyContinuation
    const HC = HomotopyContinuation
    using DocStringExtensions
    using ..HarmonicBalance

    include("HC_wrapper/homotopy_interface.jl")
end
