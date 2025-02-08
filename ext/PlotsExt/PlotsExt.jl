module PlotsExt

using DocStringExtensions
using Plots
using HarmonicBalance:
    HarmonicBalance,
    Result,
    branch_count,
    transform_solutions,
    _apply_mask,
    _get_mask,
    _realify,
    get_class,
    get_variable_solutions,
    phase_diagram,
    swept_parameters
using Symbolics: Num
using LinearAlgebra: LinearAlgebra, eigvals, eigvecs

using Plots: heatmap, theme_palette, scatter, RGB, cgrad
using Symbolics.Latexify: Latexify, latexify, @L_str
using Symbolics.Latexify.LaTeXStrings: LaTeXStrings
using SciMLBase: SciMLBase

const _set_Plots_default = Dict{Symbol,Any}([
    :fontfamily => "computer modern",
    :titlefont => "computer modern",
    :tickfont => "computer modern",
    :linewidth => 2,
    :legend_position => :outerright,
])

include("steady_states.jl")
include("linear_response.jl")
include("time_evolution.jl")

end #module
