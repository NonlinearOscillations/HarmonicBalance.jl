"""
$(TYPEDEF)

Stores the steady states of a HarmonicEquation.

# Fields
$(TYPEDFIELDS)

"""
struct Result{D,SolType<:Number,ParType<:Number,F<:JacobianFunction(SolType)}
    "The variable values of steady-state solutions."
    solutions::Array{Vector{Vector{SolType}},D}
    "Values of all parameters for all solutions."
    swept_parameters::OrderedDict{Num,Vector{ParType}}
    "The parameters fixed throughout the solutions."
    fixed_parameters::OrderedDict{Num,ParType}
    "The `Problem` used to generate this."
    problem::Problem{ParType,F}
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

    return Result{dim,soltype,partype,typeof(jacobian)}(
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

dimension(res::Result) = length(size(res.solutions)) # give solution dimensionality

"""
    phase_diagram(res::Result{D}; class="physical", not_class=[]) where {D}

Calculate the phase diagram from a `Result` object by summing over the number of
states at each swept parameters.

# Keyword arguments
Class selection done by passing `String` or `Vector{String}` as kwarg:

    class::String       :   only count solutions in this class ("all" --> plot everything)
    not_class::String   :   do not count solutions in this class

# Returns
- Array{Int64,D}: Sum of states after applying the specified class masks
"""
function phase_diagram(res::Result; class="physical", not_class=[])
    return Z = sum.(_get_mask(res, class, not_class))
end

function swept_parameters(res::Result{D}) where {D}
    X = collect(values(res.swept_parameters))
    return D == 1 ? X[1] : X
end

swept_parameter(res::Result, x::Num) = res.swept_parameters[x]
function swept_parameter(res::Result, x::String)
    is_swept_parameter(res, x)
    return res.swept_parameters[HarmonicBalance._parse_expression(x)]
end

"""
    attractors(res::Result{D}; class="stable", not_class=[]) where D

Extract attractors from a [`Result`](@ref) object. Returns an array of dictionaries, where
each dictionary maps branch identifer to the attractor. The attractors are filtered by their
corresponding class.

# Keyword arguments
Class selection done by passing `String` or `Vector{String}` as kwarg:

    class::String       :   only count solutions in this class ("all" --> plot everything)
    not_class::String   :   do not count solutions in this class


# Returns
`Array{Dict,D}`: Vector of dictionaries mapping branch indices to points satisfying
  the stability criteria at each parameter value
"""
function attractors(res::Result{D}; class="stable", not_class=[]) where {D}
    branches = 1:branch_count(res)
    Y = _get_mask(res, class, not_class)

    return map(enumerate(Y)) do (idx, bools)
        Dict(i => get_branch(res, i)[idx] for (i, bool) in pairs(bools) if bool)
    end # map
end

function is_variable(res::Result, x::String)
    vars = res.problem.variables
    x_index = findfirst(sym -> string(sym) == x, vars)
    isnothing(x_index) && error("The variable $x is not a defined variable.")
end

function is_swept_parameter(res::Result, z::String)
    # compare strings because type Num cannot be compared
    swept_pars = res.swept_parameters.keys
    z_index = findfirst(sym -> string(sym) == z, swept_pars)
    isnothing(z_index) && error("The variable $z was not swept over.")
    return true
end
is_swept_parameter(res::Result, z::Num) = is_swept_parameter(res::Result, string(z))
