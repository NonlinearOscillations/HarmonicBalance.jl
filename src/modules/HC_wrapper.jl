module HC_wrapper

using DocStringExtensions
using Symbolics: Num, @variables
using Symbolics.SymbolicUtils: isterm

using HarmonicBalance:
    HarmonicBalance,
    HarmonicEquation,
    _remove_brackets,
    expand_derivatives,
    var_name,
    get_variables,
    Problem
using HomotopyContinuation
using HomotopyContinuation: Variable, System

include("HC_wrapper/homotopy_interface.jl")

export Problem

end
