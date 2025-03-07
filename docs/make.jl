CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using HarmonicBalance

using Documenter
using DocumenterVitepress
using DocumenterCitations

# extensions
using ModelingToolkit
using OrdinaryDiffEqTsit5
using SteadyStateDiffEq

TimeEvolution = Base.get_extension(HarmonicBalance, :TimeEvolution)
ModelingToolkitExt = Base.get_extension(HarmonicBalance, :ModelingToolkitExt)
SteadyStateDiffEqExt = Base.get_extension(HarmonicBalance, :SteadyStateDiffEqExt)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric,  # default
)

using Plots
PlotsExt = Base.get_extension(HarmonicBalance, :PlotsExt)
default(; fmt=:png)
# Gotta set this environment variable when using the GR run-time on CI machines.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

include("make_md_examples.jl")

include("pages.jl")

makedocs(;
    sitename="HarmonicBalance.jl",
    authors="Quest group",
    modules=[
        HarmonicBalance,
        TimeEvolution,
        ModelingToolkitExt,
        SteadyStateDiffEqExt,
        HarmonicBalance.LinearResponse,
        PlotsExt,
    ],
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/QuantumEngineeredSystems/HarmonicBalance.jl",
        devbranch="master",
        devurl="dev",
    ),
    checkdocs=:exports,
    pages=pages,
    source="src",
    build="build",
    draft=!CI,
    warnonly=if CI
        [:linkcheck, :cross_references]
    else
        [:linkcheck, :cross_references, :missing_docs, :docs_block]
    end,
    doctest=false,  # We test it in the CI, no need to run it here
    plugins=[bib],
)

if CI
    deploydocs(;
        repo="github.com/QuantumEngineeredSystems/HarmonicBalance.jl",
        devbranch="master",
        target="build",
        branch="gh-pages",
        push_preview=true,
    )
end
