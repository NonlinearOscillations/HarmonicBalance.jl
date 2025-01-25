abstract type Problem end

"""
$(TYPEDEF)

Holds a set of algebraic equations describing the steady state of a system.

# Fields
$(TYPEDFIELDS)

#  Constructors
```julia
HomotopyContinuationProblem(
    eom::HarmonicEquation,
    swept::AbstractDict,
    fixed::AbstractDict;
    compute_Jacobian::Bool=true,
)
```
"""
mutable struct HomotopyContinuationProblem{
    Jac<:JacobianFunction(ComplexF64), # HC.jl only supports Float64
    ParType<:Number,
} <: Problem
    "The harmonic variables to be solved for."
    variables::Vector{Num}
    "All symbols which are not the harmonic variables."
    parameters::Vector{Num}
    "The swept parameters in the homotopy."
    swept_parameters::OrderedDict{Num,Vector{ParType}}
    "The fixed parameters in the homotopy."
    fixed_parameters::OrderedDict{Num,ParType}
    "The input object for HomotopyContinuation.jl solver methods."
    system::HC.System
    """
    The Jacobian matrix (possibly symbolic or compiled function).
    If `Matrix{Nan}` and implicit function is compiled when a `Result` is created.
    """
    jacobian::Jac
    "The HarmonicEquation object used to generate this `Problem`."
    eom::HarmonicEquation

    function HomotopyContinuationProblem(
        variables, parameters, swept, fixed::OrderedDict{K,V}, system, jacobian
    ) where {K,V}
        return new{typeof(jacobian),V}(
            variables, parameters, swept, fixed, system, jacobian
        )
    end # incomplete initialization for user-defined symbolic systems
    function HomotopyContinuationProblem(
        variables, parameters, swept, fixed::OrderedDict{K,V}, system, jacobian, eom
    ) where {K,V}
        return new{typeof(jacobian),V}(
            variables, parameters, swept, fixed, system, jacobian, eom
        )
    end
end

"Constructor for the type `Problem` (to be solved by HomotopyContinuation)
from a `HarmonicEquation`."
function HarmonicBalance.HomotopyContinuationProblem(
    eom::HarmonicEquation, swept::AbstractDict, fixed::AbstractDict
)
    S = HomotopyContinuation.System(eom)
    vars_new = declare_variables(eom)

    swept, fixed = promote_types(swept, fixed)
    check_fixed_and_sweep(eom, swept, fixed)

    jac = _compile_Jacobian(eom, ComplexF64, swept, fixed)
    # ^ HC.jl only supports Float64 (https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl/issues/604)
    return HomotopyContinuationProblem(vars_new, eom.parameters, swept, fixed, S, jac, eom)
end

Symbolics.get_variables(p::Problem)::Vector{Num} = get_variables(p.eom)

function Base.show(io::IO, p::HomotopyContinuationProblem)
    println(io, length(p.system.expressions), " algebraic equations for steady states")
    println(io, "Variables: ", join(string.(p.variables), ", "))
    println(io, "Parameters: ", join(string.(p.parameters), ", "))
    return println(io, "Symbolic Jacobian: ", !(p.jacobian == false))
end

# assume this order of variables in all compiled function (transform_solutions, Jacobians)
function _free_symbols(p::Problem)::Vector{Num}
    return cat(p.variables, collect(keys(p.swept_parameters)); dims=1)
end

function declare_variables(p::Problem)
    return declare_variable.(string.(cat(p.parameters, p.variables; dims=1)))
end

function check_fixed_and_sweep(
    eom::Union{HomotopyContinuationProblem,HarmonicEquation}, sweeps, fixed_parameters
)
    variable_names = var_name.([keys(fixed_parameters)..., keys(sweeps)...])
    if any([var_name(var) ∈ variable_names for var in get_variables(eom)])
        error("Cannot fix one of the variables!")
    end
    # sweeping takes precedence over fixed_parameters
    unique_fixed = filter_duplicate_parameters(sweeps, fixed_parameters)

    # fixed order of parameters
    all_keys = cat(collect(keys(sweeps)), collect(keys(fixed_parameters)); dims=1)
    # the order of parameters we have now does not correspond to that in prob!
    # get the order from prob and construct a permutation to rearrange our parameters
    error = ArgumentError("Some input parameters are missing or appear more than once!")
    permutation = try
        # find the matching position of each parameter
        p = [findall(x -> isequal(x, par), all_keys) for par in eom.parameters]
        all((length.(p)) .== 1) || throw(error) # some parameter exists more than twice!
        p = getindex.(p, 1) # all exist once -> flatten
        # ∨ parameters sorted wrong!
        isequal(all_keys[getindex.(p, 1)], eom.parameters) || throw(error)
        # ∨ extra parameters present!
        #isempty(setdiff(all_keys, prob.parameters)) || throw(error)
        p
    catch
        throw(error)
    end
    return unique_fixed, permutation
end
