module LinearResponse

using Printf: Printf, @printf
using DocStringExtensions
using ProgressMeter: ProgressMeter, Progress, next!

using Symbolics: Symbolics, Num, unwrap
using LinearAlgebra: norm, eigen

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

export get_Jacobian,
    show, get_jacobian_response, get_linear_response, get_rotframe_jacobian_response

end
