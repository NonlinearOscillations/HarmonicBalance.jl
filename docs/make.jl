CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using HarmonicBalance
using Documenter
using DocumenterVitepress
using DocumenterCitations

# extentions
using ModelingToolkit
using OrdinaryDiffEq
using SteadyStateDiffEq

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric,  # default
)

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
