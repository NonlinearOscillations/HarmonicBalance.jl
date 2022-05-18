import HomotopyContinuation: Variable
import HarmonicBalance: Problem; export Problem
export Num_to_Variable

"Conversion from Symbolics.jl types to HomotopyContinuation types."
Variable(var::Num) = var.val isa SymbolicUtils.Term ? Variable(string(var.val.f)) : Variable(string(var_name(var)))


"Converts a Num into Variable in the active namespace."
function Num_to_Variable(x::Num)
    var = Variable(x)
    s = Symbol(string(var))
    @eval (($s) = ($var))
end

"Converts a Num dictionary into a Variable dictionary."
Num_to_Variable(dict::Dict{Num, ComplexF64}) = Dict{Variable, ComplexF64}([[Variable(key), dict[key]] for key in keys(dict)]) # for the parameter assignments


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
function Problem(eom::HarmonicEquation; Jacobian=true)

    S = System(eom)
    # use the rearranged system for the proper definition of the Jacobian
    # this possibly has variables in the denominator and cannot be used for solving
    if Jacobian == true || Jacobian == "explicit"
        J = HarmonicBalance.get_Jacobian(eom)
    elseif Jacobian == "implicit"
        # compute the Jacobian implicitly
        J = HarmonicBalance.LinearResponse.get_implicit_Jacobian(eom)
    else
        J = Jacobian
    end
    vars_orig  = get_variables(eom)
    vars_new = declare_variable.(HarmonicBalance.var_name.(vars_orig))
    return Problem(vars_new, eom.parameters, S, J,eom)   
end


"A constructor for Problem from explicitly entered equations, variables and parameters."
function Problem(equations::Vector{Num},variables::Vector{Num},parameters::Vector{Num})
    conv_vars = Num_to_Variable.(variables)
    conv_para = Num_to_Variable.(parameters)
    
    eqs_HC=[Expression(eval(symbol)) for symbol in [Meta.parse(s) for s in [string(eq) for eq in equations]]] #note in polar coordinates there could be imaginary factors, requiring the extra replacement "I"=>"1im"
    system = HomotopyContinuation.System(eqs_HC, variables = conv_vars, parameters = conv_para)
    J = HarmonicBalance.get_Jacobian(equations,variables) #all derivatives are assumed to be in the left hand side;
    return Problem(variables,parameters,system,J)
end


function System(eom::HarmonicEquation)
    eqs = expand_derivatives.(_remove_brackets(eom))
    conv_vars = Num_to_Variable.(get_variables(eom))
    conv_para = Num_to_Variable.(eom.parameters)
    S = HomotopyContinuation.System(parse_equations(eqs),variables=conv_vars,parameters=conv_para)
end


###
# DEPRECATED
###

#=

"Parse symbolic expressions as the Expression type in HomotopyContinuation."
function parse_equations(eqs::HarmonicEquation)
    equation_strings = string.(HarmonicBalance._equations_without_brackets(eqs))
    println(equation_strings)
    parsed_strings = [Meta.parse(s) for s in equation_strings]
    return [Expression(eval(symbol)) for symbol in parsed_strings]
end

"""
$(TYPEDSIGNATURES)
Return the Jacobian of the solution `s` (specifies all variables and parameters)
of `problem`.
"""
function Jacobian(s::Dict, problem::Problem; im_tol=im_tol)
    #svar = Num_to_Variable(s)
    #isnan(prod(values(svar))) && error("trying to find the Jacobian of an invalid solution")
    #!HarmonicBalance.is_physical(s, problem, im_tol=im_tol) && error("unphysical solutions have no Jacobian")
    return substitute_all(problem.jacobian, s)
end

function Jacobian(s::Dict, res::HarmonicBalance.Result; im_tol=im_tol)
    J = res.jacobian
    vals = [s[var] for var in cat(res.problem.variables, collect(keys(res.swept_parameters)), dims=1)]
    return real.([Base.invokelatest(el, vals) for el in J])
end

function construct_Jacobian(eom::HarmonicEquation)
    parameters = declare_variable.(eom.parameters)
    time = declare_variable(get_independent_variables(eom)[1])
    variables = declare_variable.(HarmonicBalance.var_name.(get_variables(eom)))
    dummy_system = deepcopy(eom)
    eqs = HarmonicBalance.equation_without_brackets(eom)
    #dummy_system.equations = [eval(Meta.parse(string(eq))) for eq in eqs]
    #Symbolics.solve_for(dummy_system.equations .~ dummy_system.left_hand_side, d.(get_variables(eom), time))
    #S_rearranged = HomotopyContinuation.System(equations_hc(rearranged),variables=conv_vars,parameters=conv_para) 
end
=#

#=
"Converts a Variable into Sym in the HomotopySolving namespace"
function Variable_to_Num(x::Variable)
    #s = Symbol(x.name)
    s = Symbol(string(x))
    var = @variables($s)
    var = var[1]
    display(typeof(var))
    s = Symbol(string(x))
    @eval (($s) = ($var))   
end

"Performs multiple replacements on a string. REPLACEMENTS ARE PERFORMED IN SUCCESSION, EACH ONLY ONCE"
function replace_string_successive(old_string::String, rules::Vector{Pair{String, String}})
    new_string = old_string
    for rule in rules
        new_string = Base.replace(new_string, rule)
    end
    new_string
end
=#


