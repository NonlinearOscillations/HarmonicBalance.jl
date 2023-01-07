module TimeEvolution

    using ..HarmonicBalance
    using Symbolics
    using OrdinaryDiffEq
    using DSP
    using FFTW
    using Peaks
    using DocStringExtensions

    include("TimeEvolution/types.jl")
    include("TimeEvolution/ODEProblem.jl")
    include("TimeEvolution/FFT_analysis.jl")
    include("TimeEvolution/sweeps.jl")

end
