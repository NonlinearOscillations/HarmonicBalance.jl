struct JacobianFunction{Ret,Args<:Tuple}
    J::FunctionWrapper{Ret,Args}
end

(cb::JacobianFunction)(v) = cb.J(v)


""" Compile the Jacobian from `prob`, inserting `fixed_parameters`.
    Returns a function that takes a dictionary of variables and `swept_parameters` to give the Jacobian."""
function _compile_Jacobian(
    prob::Problem, swept_parameters::OrderedDict, fixed_parameters::OrderedDict
)
    if prob.jacobian isa Matrix
        compiled_J = compile_matrix(
            prob.jacobian, _free_symbols(prob, swept_parameters); rules=fixed_parameters
        )
    elseif prob.jacobian == "implicit"
        compiled_J = LinearResponse.get_implicit_Jacobian(
            prob, swept_parameters, fixed_parameters
        ) # leave implicit Jacobian as is
    else
        return prob.jacobian
    end
    return compiled_J
end


"""
Take a matrix containing symbolic variables `variables` and keys of `fixed_parameters`.
Substitute the values according to `fixed_parameters` and compile into a function that takes numerical arguments
    in the order set in `variables`.
"""
function compile_matrix(mat, variables; rules=Dict(), postproc=x -> x)
    J = substitute_all.(mat, Ref(rules))
    matrix = Symbolics.build_function(J, variables)
    matrix = eval(matrix[1]) # compiled allocating function, see Symbolics manual
    m(vals::Vector) = postproc(matrix(vals))
    m(s::OrderedDict) = m([s[var] for var in variables]) # for the UI
    return m
end
