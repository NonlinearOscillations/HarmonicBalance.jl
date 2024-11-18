
"""
$(TYPEDEF)

Holds a set of algebraic equations describing the steady state of a system.

# Fields
$(TYPEDFIELDS)

#  Constructors
```julia
Problem(eom::HarmonicEquation; Jacobian=true) # find and store the symbolic Jacobian
Problem(eom::HarmonicEquation; Jacobian="implicit") # ignore the Jacobian for now, compute implicitly later
Problem(eom::HarmonicEquation; Jacobian=J) # use J as the Jacobian (a function that takes a Dict)
Problem(eom::HarmonicEquation; Jacobian=false) # ignore the Jacobian
```
"""
struct Problem
    "The harmonic variables to be solved for."
    variables::Vector{Num}
    "All symbols which are not the harmonic variables."
    parameters::Vector{Num}
    "The input object for HomotopyContinuation.jl solver methods."
    system::HC.System
    "The Jacobian matrix (possibly symbolic).
    If `false`, the Jacobian is ignored (may be calculated implicitly after solving)."
    jacobian
    "The HarmonicEquation object used to generate this `Problem`."
    eom::HarmonicEquation

    function Problem(variables, parameters, system, jacobian)
        return new(variables, parameters, system, jacobian)
    end #incomplete initialization for user-defined symbolic systems
    function Problem(variables, parameters, system, jacobian, eom)
        return new(variables, parameters, system, jacobian, eom)
    end
end

Symbolics.get_variables(p::Problem)::Vector{Num} = get_variables(p.eom)

function Base.show(io::IO, p::Problem)
    println(io, length(p.system.expressions), " algebraic equations for steady states")
    println(io, "Variables: ", join(string.(p.variables), ", "))
    println(io, "Parameters: ", join(string.(p.parameters), ", "))
    return println(io, "Symbolic Jacobian: ", !(p.jacobian == false))
end

# assume this order of variables in all compiled function (transform_solutions, Jacobians)
function _free_symbols(p::Problem, varied)
    return cat(p.variables, collect(keys(OrderedDict(varied))); dims=1)
end

""" Compile the Jacobian from `prob`, inserting `fixed_parameters`.
    Returns a function that takes a dictionary of variables and `swept_parameters` to give the Jacobian."""
function _compile_Jacobian(
    prob::Problem, swept_parameters::ParameterRange, fixed_parameters::ParameterList
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
