module PlotsExt

using DocStringExtensions
using Plots
using HarmonicBalance
using HarmonicBalance:
    HarmonicBalance,
    Result,
    branch_count,
    transform_solutions,
    _realify,
    get_class,
    get_variable_solutions,
    phase_diagram,
    is_variable,
    is_swept_parameter,
    _get_mask

using Symbolics: Num
using LinearAlgebra: LinearAlgebra, eigvals, eigvecs

using HarmonicBalance.LinearResponse:
    get_jacobian_response, get_linear_response, get_rotframe_jacobian_response

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

function get_labels(res::Result{D}) where D
    if D == 1
        return latexify(string(first(keys(res.swept_parameters))))
    elseif D ==2
        return latexify.(string.(keys(res.swept_parameters)))
    else
        error("Getting the labels is only supported for 1D and 2D results.")
    end
end

include("steady_states.jl")
include("linear_response.jl")
include("time_evolution.jl")

end #module
