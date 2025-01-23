module LinearResponse

using Printf: Printf, @printf
using DocStringExtensions
using ProgressMeter: ProgressMeter, Progress, next!

using Plots: heatmap, theme_palette, scatter, RGB, cgrad
using Latexify: Latexify, latexify, @L_str
using Latexify.LaTeXStrings: LaTeXStrings

using Symbolics: Symbolics, Num, unwrap
using LinearAlgebra: norm, eigvals, eigen, eigvecs

using HarmonicBalance
using HarmonicBalance:
    Result,
    HarmonicVariable,
    DifferentialEquation,
    StateDict,
    get_variables,
    get_independent_variables,
    get_variable_solutions

using HarmonicBalance:
    var_name,
    d,
    substitute_all,
    _free_symbols,
    _get_mask,
    _set_Plots_default,
    dim,
    _get_mask,
    harmonic_ansatz,
    slow_flow,
    fourier_transform,
    declare_variable

using ..HC_wrapper

include("types.jl")
include("utils.jl")
include("Lorentzian_spectrum.jl")
include("response.jl")
include("plotting.jl")

export get_Jacobian, show
export plot_linear_response, plot_rotframe_jacobian_response, plot_eigenvalues

end
