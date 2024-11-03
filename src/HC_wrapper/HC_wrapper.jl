module HC_wrapper

using DocStringExtensions
using Symbolics: Num, @variables, expand_derivatives, get_variables
using Symbolics.SymbolicUtils: isterm
using LinearAlgebra: LinearAlgebra

using HarmonicBalance:
    HarmonicBalance, HarmonicEquation, _remove_brackets, var_name, Problem

using HomotopyContinuation
using HomotopyContinuation: Variable, System

include("homotopy_interface.jl")

export Problem

end
