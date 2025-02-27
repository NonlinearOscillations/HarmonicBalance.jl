module HC_wrapper

using DocStringExtensions
using Symbolics: Num, @variables, expand_derivatives, get_variables
using Symbolics.SymbolicUtils: isterm
using LinearAlgebra: LinearAlgebra

using HarmonicBalance:
    HarmonicBalance, HarmonicEquation, _remove_brackets, var_name, Problem

using HomotopyContinuation
using HomotopyContinuation: Variable, System

"Conversion from Symbolics.jl types to HomotopyContinuation types."
HomotopyContinuation.Variable(var::Num) =
    isterm(var.val) ? Variable(string(var.val.f)) : Variable(string(var_name(var)))

"Converts a Num into Variable in the active namespace."
function Num_to_Variable(x::Num)
    var = Variable(x)
    s = Symbol(string(var))
    @eval (($s) = ($var))
end

"Converts a Num dictionary into a Variable dictionary."
Num_to_Variable(dict::Dict{Num,T}) where {T<:Number} =
    Dict{Variable,T}([[Variable(key), dict[key]] for key in keys(dict)]) # for the parameter assignments

"Parse symbolic expressions as the Expression type in HomotopyContinuation."
function parse_equations(eqs::Vector{Num})
    parsed_strings = [Meta.parse(s) for s in string.(eqs)]
    return [Expression(eval(symbol)) for symbol in parsed_strings]
end

function HomotopyContinuation.System(eom::HarmonicEquation)
    eqs = expand_derivatives.(_remove_brackets(eom))
    vars = get_variables(eom)
    pars = eom.parameters
    return  System(eqs, vars, pars)
end
function HomotopyContinuation.System(eqs::Vector{Num}, vars::Vector{Num}, pars::Vector{Num})
    conv_vars = Num_to_Variable.(vars)
    conv_para = Num_to_Variable.(pars)
    return S = HomotopyContinuation.System(
        parse_equations(eqs); variables=conv_vars, parameters=conv_para
    )
end

export Problem

end
