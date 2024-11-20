"""
$(TYPEDEF)

Stores the steady states of a HarmonicEquation.

# Fields
$(TYPEDFIELDS)

"""
struct Result{SolType<:Number,ParType<:Number,D,F<:JacobianFunction(SolType),F1}
    "The variable values of steady-state solutions."
    solutions::Array{Vector{Vector{SolType}},D}
    "Values of all parameters for all solutions."
    swept_parameters::OrderedDict{Num,Vector{ParType}}
    "The parameters fixed throughout the solutions."
    fixed_parameters::OrderedDict{Num,ParType}
    "The `Problem` used to generate this."
    problem::Problem{F1}
    """
    Maps strings such as \"stable\", \"physical\" etc to arrays of values,
    classifying the solutions (see method `classify_solutions!`).
    """
    classes::Dict{String,Array{BitVector,D}}
    """
    Create binary classification of the solutions, such that each solution point receives
    an identifier based on its permutation of stable branches (allows to distinguish between
    different phases, which may have the same number of stable solutions). It works by
    converting each bitstring `[is_stable(solution_1), is_stable(solution_2), ...,]` into
    unique labels.
    """
    binary_labels::Array{Int64,D}
    """
    The Jacobian function with `fixed_parameters` already substituted. Accepts a vector
    specifying the solution. If problem.jacobian is a symbolic matrix, this holds a compiled
    function.
    """
    jacobian::F
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
    soltype = solution_type(solutions)
    partype = parameter_type(swept_parameters, fixed_parameters)
    dim = ndims(solutions)

    return Result{soltype,partype,dim,typeof(jacobian),typeof(problem.jacobian)}(
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
    println(io, "   of which real:    ", sum(push!(any.(get_class(r, "physical")), false)))
    println(io, "   of which stable:  ", sum(push!(any.(get_class(r, "stable")), false)))
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
