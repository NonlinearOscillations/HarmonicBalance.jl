"""
$(TYPEDEF)

Stores the steady states of a HarmonicEquation.

# Fields
$(TYPEDFIELDS)

"""
struct Result{SolutionType<:Number,ParameterType<:Number,VarLength,Dimension}
    "The variable values of steady-state solutions."
    solutions::Array{Vector{Vector{SolutionType}},Dimension}
    "Values of all parameters for all solutions."
    swept_parameters::OrderedDict{Num,Vector{ParameterType}}
    "The parameters fixed throughout the solutions."
    fixed_parameters::OrderedDict{Num,ParameterType}
    "The `Problem` used to generate this."
    problem::Problem
    "Maps strings such as \"stable\", \"physical\" etc to arrays of values, classifying the solutions (see method `classify_solutions!`)."
    classes::Dict{String,Array{BitVector,Dimension}}
    "Create binary classification of the solutions, such that each solution point receives an identifier
    based on its permutation of stable branches (allows to distinguish between different phases,
    which may have the same number of stable solutions). It works by converting each bitstring
    `[is_stable(solution_1), is_stable(solution_2), ...,]` into unique labels."
    binary_labels::Array{Int64,Dimension}
    "The Jacobian with `fixed_parameters` already substituted. Accepts a dictionary specifying the solution.
    If problem.jacobian is a symbolic matrix, this holds a compiled function.
    If problem.jacobian was `false`, this holds a function that rearranges the equations to find J
    only after numerical values are inserted (preferable in cases where the symbolic J would be very large)."
    jacobian::JacobianFunction{Matrix{ParameterType},NTuple{VarLength,SolutionType}}
    "Seed used for the solver"
    seed::UInt32
end

function Result(
    solutions,
    swept_parameters,
    fixed_parameters,
    problem,
    classes,
    binary_labels,
    jacobian,
    seed,
)
    soltype = eltype(eltype(eltype(solutions)))
    partype = eltype(eltype(swept_parameters).parameters[2])
    dim = ndims(solutions)
    varlength = length(swept_parameters)+length(get_variables(problem))

    return Result{soltype,partype,varlength,dim}(
        solutions,
        swept_parameters,
        fixed_parameters,
        problem,
        classes,
        binary_labels,
        jacobian,
        seed,
    )
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
