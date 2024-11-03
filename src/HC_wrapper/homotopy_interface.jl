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
Num_to_Variable(dict::Dict{Num,ComplexF64}) =
    Dict{Variable,ComplexF64}([[Variable(key), dict[key]] for key in keys(dict)]) # for the parameter assignments

"Parse symbolic expressions as the Expression type in HomotopyContinuation."
function parse_equations(eqs::Vector{Num})
    parsed_strings = [Meta.parse(s) for s in string.(eqs)]
    return [Expression(eval(symbol)) for symbol in parsed_strings]
end

"Declare a new variable in the the current namespace."
function declare_variable(name::String)
    var_sym = Symbol(name)
    @eval($(var_sym) = first(@variables $var_sym))
    return eval(var_sym)
end

declare_variable(x::Num) = declare_variable(string(x))

"Constructor for the type `Problem` (to be solved by HomotopyContinuation)
from a `HarmonicEquation`."
function HarmonicBalance.Problem(eom::HarmonicEquation; Jacobian=true)
    S = System(eom)
    # use the rearranged system for the proper definition of the Jacobian
    # this possibly has variables in the denominator and cannot be used for solving
    if Jacobian == true || Jacobian == "explicit"
        J = HarmonicBalance.get_Jacobian(eom)
    elseif Jacobian == "false" || Jacobian == false
        dummy_J(arg) = LinearAlgebra.I(1)
        J = dummy_J
    else
        J = Jacobian
    end
    vars_orig = get_variables(eom)
    vars_new = declare_variable.(var_name.(vars_orig))
    return Problem(vars_new, eom.parameters, S, J, eom)
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
end # TODO is this funciton still needed/used?

function System(eom::HarmonicEquation)
    eqs = expand_derivatives.(_remove_brackets(eom))
    conv_vars = Num_to_Variable.(get_variables(eom))
    conv_para = Num_to_Variable.(eom.parameters)
    return S = HomotopyContinuation.System(
        parse_equations(eqs); variables=conv_vars, parameters=conv_para
    )
end
