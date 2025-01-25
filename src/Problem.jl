
"""
$(TYPEDEF)

Holds a set of algebraic equations describing the steady state of a system.

# Fields
$(TYPEDFIELDS)

#  Constructors
```julia
Problem(eom::HarmonicEquation)
```
"""
mutable struct Problem{Jac}
    "The harmonic variables to be solved for."
    variables::Vector{Num}
    "All symbols which are not the harmonic variables."
    parameters::Vector{Num}
    "The input object for HomotopyContinuation.jl solver methods."
    system::HC.System
    """
    The Jacobian matrix (possibly symbolic or compiled function).
    If `Matrix{Nan}` and implicit function is compiled when a `Result` is created.
    """
    jacobian::Jac
    "The HarmonicEquation object used to generate this `Problem`."
    eom::HarmonicEquation

    function Problem(variables, parameters, system, jacobian)
        return new{typeof(jacobian)}(variables, parameters, system, jacobian)
    end # incomplete initialization for user-defined symbolic systems
    function Problem(variables, parameters, system, jacobian, eom)
        return new{typeof(jacobian)}(variables, parameters, system, jacobian, eom)
    end
end

"Constructor for the type `Problem` (to be solved by HomotopyContinuation)
from a `HarmonicEquation`."
function HarmonicBalance.Problem(eom::HarmonicEquation; compute_Jacobian::Bool=true)
    S = HomotopyContinuation.System(eom)
    vars_orig = get_variables(eom)
    vars_new = declare_variable.(var_name.(vars_orig))

    jac = compute_Jacobian ? eom.Jacobian : dummy_symbolic_Jacobian(length(eom.variables))
    return Problem(vars_new, eom.parameters, S, jac, eom)
end

Symbolics.get_variables(p::Problem)::Vector{Num} = get_variables(p.eom)

function Base.show(io::IO, p::Problem)
    println(io, length(p.system.expressions), " algebraic equations for steady states")
    println(io, "Variables: ", join(string.(p.variables), ", "))
    println(io, "Parameters: ", join(string.(p.parameters), ", "))
    return println(io, "Symbolic Jacobian: ", !(p.jacobian == false))
end

# assume this order of variables in all compiled function (transform_solutions, Jacobians)
function _free_symbols(p::Problem, varied::OrderedDict)
    return cat(p.variables, collect(keys(varied)); dims=1)
end
