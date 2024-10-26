module LinearResponse

using Printf: Printf, @printf
using DocStringExtensions
using ProgressMeter: ProgressMeter, Progress, next!

using Plots: heatmap, theme_palette, scatter, RGB, cgrad
using Latexify: Latexify, latexify, @L_str
using Latexify.LaTeXStrings: LaTeXStrings

using Symbolics: Symbolics, Num, Equation, unwrap
using LinearAlgebra: norm, eigvals, eigen, eigvecs
using OrderedCollections: OrderedDict

using HarmonicBalance
using HarmonicBalance:
    var_name,
    rearrange_standard,
    _remove_brackets,
    expand_derivatives,
    substitute_all,
    _free_symbols,
    _get_mask,
    compile_matrix,
    _set_Plots_default,
    dim,
    _get_mask,
    harmonic_ansatz,
    slow_flow,
    fourier_transform,
    declare_variable,
    is_rearranged
using ..HC_wrapper

include("types.jl")
include("utils.jl")
include("jacobians.jl")
include("Lorentzian_spectrum.jl")
include("response.jl")
include("plotting.jl")

export get_Jacobian, get_response, show
export plot_linear_response, plot_rotframe_jacobian_response, plot_eigenvalues

end
