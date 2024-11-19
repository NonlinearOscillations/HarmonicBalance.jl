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

"A constructor for Problem from explicitly entered equations, variables and parameters."
function HarmonicBalance.Problem(
    equations::Vector{Num}, variables::Vector{Num}, parameters::Vector{Num}
)
    conv_vars = Num_to_Variable.(variables)
    conv_para = Num_to_Variable.(parameters)

    eqs_HC = [
        Expression(eval(symbol)) for
        symbol in [Meta.parse(s) for s in [string(eq) for eq in equations]]
    ] #note in polar coordinates there could be imaginary factors, requiring the extra replacement "I"=>"1im"
    system = HomotopyContinuation.System(eqs_HC; variables=conv_vars, parameters=conv_para)
    J = HarmonicBalance.get_Jacobian(equations, variables) #all derivatives are assumed to be in the left hand side;
    return Problem(variables, parameters, system, J)
end # TODO is this function still needed/used?

function System(eom::HarmonicEquation)
    eqs = expand_derivatives.(_remove_brackets(eom))
    conv_vars = Num_to_Variable.(get_variables(eom))
    conv_para = Num_to_Variable.(eom.parameters)
    return S = HomotopyContinuation.System(
        parse_equations(eqs); variables=conv_vars, parameters=conv_para
    )
end

export Problem

end
