push!(LOAD_PATH, "../src/")

using Documenter
using HarmonicBalance
using ModelingToolkit
using OrdinaryDiffEq
using SteadyStateDiffEq

include("pages.jl")

makedocs(;
    sitename="HarmonicBalance.jl",
    authors="Nonlinear Oscillations Group",
    modules=[
        HarmonicBalance,
        Base.get_extension(HarmonicBalance, :TimeEvolution),
        Base.get_extension(HarmonicBalance, :ModelingToolkitExt),
        Base.get_extension(HarmonicBalance, :SteadyStateDiffEqExt),
    ],
    warnonly=true,
    format=Documenter.HTML(;
        mathengine=MathJax2(),
        canonical="https://nonlinearoscillations.github.io/HarmonicBalance.jl/stable/",
        assets=["assets/favicon.ico", "assets/docs.css"],
        # size_threshold = nothing
    ),
    pages=pages,
)

deploydocs(;
    repo="github.com/NonlinearOscillations/HarmonicBalance.jl.git", push_preview=false
)
