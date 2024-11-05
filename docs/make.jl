CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using HarmonicBalance
using Documenter
using DocumenterVitepress
using DocumenterCitations

# extensions
using ModelingToolkit
using OrdinaryDiffEqTsit5
using SteadyStateDiffEq

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric,  # default
)

using Plots
default(; fmt=:png)
# Gotta set this environment variable when using the GR run-time on CI machines.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

include("make_md_examples.jl")

makedocs(;
    sitename="HarmonicBalance.jl",
    authors="Quest",
    modules=[
        HarmonicBalance,
        Base.get_extension(HarmonicBalance, :TimeEvolution),
        Base.get_extension(HarmonicBalance, :ModelingToolkitExt),
        Base.get_extension(HarmonicBalance, :SteadyStateDiffEqExt),
    ],
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/NonlinearOscillations/HarmonicBalance.jl",
        devbranch="master",
        devurl="dev",
    ),
    source="src",
    build="build",
    draft=false,
    warnonly=true,
    doctest=false,  # We test it in the CI, no need to run it here
    plugins=[bib],
)

if CI
    deploydocs(;
        repo="github.com/NonlinearOscillations/HarmonicBalance.jl",
        devbranch="master",
        target="build",
        branch="gh-pages",
        push_preview=true,
    )
end
