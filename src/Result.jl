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
    classes::Dict{String,Array}
    "The Jacobian with `fixed_parameters` already substituted. Accepts a dictionary specifying the solution.
    If problem.jacobian is a symbolic matrix, this holds a compiled function.
    If problem.jacobian was `false`, this holds a function that rearranges the equations to find J
    only after numerical values are inserted (preferable in cases where the symbolic J would be very large)."
    jacobian::Function
    "Seed used for the solver"
    seed::UInt32

    function Result(sol, swept, fixed, problem, classes, J, seed)
        return new(sol, swept, fixed, problem, classes, J, seed)
    end
    Result(sol, swept, fixed, problem, classes) = new(sol, swept, fixed, problem, classes)
    Result(sol, swept, fixed, problem) = new(sol, swept, fixed, problem, Dict([]))
end

Symbolics.get_variables(res::Result)::Vector{Num} = get_variables(res.problem)

function Base.show(io::IO, r::Result)
    println(io, "A steady state result for ", length(r.solutions), " parameter points")
    println(io, "\nSolution branches:   ", length(r.solutions[1]))
    println(
        io, "   of which real:    ", sum(push!(any.(classify_branch(r, "physical")), false))
    )
    println(
        io, "   of which stable:  ", sum(push!(any.(classify_branch(r, "stable")), false))
    )
    return println(io, "\nClasses: ", join(keys(r.classes), ", "))
end

# assume this order of variables in all compiled function (transform_solutions, Jacobians)
function _free_symbols(res::Result)
    return cat(res.problem.variables, collect(keys(res.swept_parameters)); dims=1)
end

# overload to use [] for indexing
Base.getindex(r::Result, idx::Int...) = get_single_solution(r, idx)
Base.size(r::Result) = size(r.solutions)

branch_count(r::Result) = length(r.solutions[1])
get_branch(r::Result, idx) = getindex.(r.solutions, idx)
