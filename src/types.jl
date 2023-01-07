using Symbolics
import HomotopyContinuation

export Problem, Result


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
mutable struct Problem
    "The harmonic variables to be solved for."
    variables::Vector{Num}
    "All symbols which are not the harmonic variables."
    parameters::Vector{Num}
    "The input object for HomotopyContinuation.jl solver methods."
    system::HomotopyContinuation.System
    "The Jacobian matrix (possibly symbolic).
    If `false`, the Jacobian is ignored (may be calculated implicitly after solving)."
    jacobian
    "The HarmonicEquation object used to generate this `Problem`."
    eom::HarmonicEquation

    Problem(variables,parameters,system,jacobian) = new(variables,parameters,system,jacobian) #incomplete initialization for user-defined symbolic systems
    Problem(variables,parameters,system,jacobian,eom) = new(variables,parameters,system,jacobian,eom)
end


function show(io::IO, p::Problem)
    println(io, length(p.system.expressions), " algebraic equations for steady states")
    println(io, "Variables: ", join(string.(p.variables), ", "))
    println(io, "Parameters: ", join(string.(p.parameters), ", "))
    println(io, "Symbolic Jacobian: ", !(p.jacobian==false))
end


"""
$(TYPEDEF)

Stores the steady states of a HarmonicEquation.

# Fields
$(TYPEDFIELDS)

"""
mutable struct Result
    "The variable values of steady-state solutions."
    solutions::Array{Vector{SteadyState}}
    "Values of all parameters for all solutions."
    swept_parameters::ParameterRange
    "The parameters fixed throughout the solutions."
    fixed_parameters::ParameterList
    "The `Problem` used to generate this."
    problem::Problem
    "Maps strings such as \"stable\", \"physical\" etc to arrays of values, classifying the solutions (see method `classify_solutions!`)."
    classes::Dict{String, Array}
    "The Jacobian with `fixed_parameters` already substituted. Accepts a dictionary specifying the solution.
    If problem.jacobian is a symbolic matrix, this holds a compiled function.
    If problem.jacobian was `false`, this holds a function that rearranges the equations to find J
    only after numerical values are inserted (preferable in cases where the symbolic J would be very large)."
    jacobian::Union{Function, Int64}

    Result(sol,swept, fixed, problem, classes, J) = new(sol, swept, fixed, problem, classes, J)
    Result(sol,swept, fixed, problem, classes) = new(sol, swept, fixed, problem, classes)
    Result(sol,swept, fixed, problem) = new(sol, swept, fixed, problem, Dict([]))
end


function show(io::IO, r::Result)
    println(io, "A steady state result for ", length(r.solutions), " parameter points")
    println(io, "\nSolution branches:   ", length(r.solutions[1]))
    println(io, "   of which real:    ", sum(any.(classify_branch(r, "physical"))))
    println(io, "   of which stable:  ", sum(any.(classify_branch(r, "stable"))))
    println(io, "\nClasses: ", join(keys(r.classes), ", "))
end


# overload to use [] for indexing
Base.getindex(r::Result, idx::Int...) = get_single_solution(r, idx)
Base.size(r::Result) = size(r.solutions)

branch_count(r::Result) = length(r.solutions[1])
get_branch(r::Result, idx) = getindex.(r.solutions, idx)
