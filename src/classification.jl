"""
$(TYPEDSIGNATURES)

Creates a solution class in `res` using the function `func` (parsed into Symbolics.jl input).
The new class is labeled with `name` and stored under `res.classes[name]`.
By default, only physical (real) solutions are classified, and `false` is returned for the rest.
To also classify complex solutions, set `physical=false`.

# Example
```julia
# solve a previously-defined problem
res = get_steady_states(problem, swept_parameters, fixed_parameters)

# classify, store in result.classes["large_amplitude"]
classify_solutions!(res, "sqrt(u1^2 + v1^2) > 1.0" , "large_amplitude")
```

"""
function classify_solutions!(
    res::Result, func::Union{String,Function}, name::String; physical=true
)
    values = classify_solutions(res, func; physical)
    return res.classes[name] = values
end

"""
$(TYPEDSIGNATURES)

Classifies solutions in `res` using the function `func`.
If `physical` is true, only physical solutions are classified.
"""
function classify_solutions(res::Result, func; physical=true)
    func = isa(func, Function) ? func : _build_substituted(func, res)
    if physical
        f_comp(soln) = _is_physical(soln) && func(real.(soln))
        return transform_solutions(res, f_comp)

    else
        return transform_solutions(res, func)
    end
end

"""
$(TYPEDSIGNATURES)

Returns an array of booleans classifying `branch` in the solutions in `res`
according to `class`.
"""
function get_class(res::Result, branch::Int64, class::String)
    return branch_values = getindex.(res.classes[class], branch)
end

"""
$(TYPEDSIGNATURES)

Returns an array of booleans classifying each branch in the solutions in `res`
according to `class`.
"""
function get_class(soln::Result, class::String)
    return [get_class(soln, b, class) for b in 1:length(first(soln.solutions))]
end

"""
$(TYPEDSIGNATURES)

Returns true if the solution `soln` of the Result `res` is physical (real number).
`im_tol` : an absolute threshold to distinguish real/complex numbers.
"""
function is_physical(soln::StateDict, res::Result)
    var_values = [getindex(soln, v) for v in res.problem.variables]
    return _is_physical(var_values)
end

_is_physical(soln; im_tol=IM_TOL) = all(x -> !isnan(x) && abs(imag(x)) < im_tol, soln)
_is_physical(res::Result) = classify_solutions(res, _is_physical)

"""
$(TYPEDSIGNATURES)

Returns true if the solution `soln` of the Result `res` is stable.
Stable solutions are real and have all Jacobian eigenvalues Re(λ) <= 0.
`im_tol` : an absolute threshold to distinguish real/complex numbers.
`rel_tol`: Re(λ) considered <=0 if real.(λ) < rel_tol*abs(λmax)
"""
function is_stable(soln::StateDict, res::Result; kwargs...)
    return _is_stable(collect(values(soln)), res.jacobian; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Returns a function that checks if a solution is stable.
"""
function _is_stable(res::Result; kwargs...)
    return _isit(soln) = _is_stable(soln, res.jacobian; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Returns true if the solution `soln` is stable.
"""
function _is_stable(soln, J; rel_tol=1e-10)
    _is_physical(soln) || return false
    λs = eigvals(real.(J(soln)))
    scale = maximum(Iterators.map(abs, λs))
    return all(x -> real(x) < rel_tol * scale, λs)
end

"""
$(TYPEDSIGNATURES)

Returns true if the solution `soln` of the problem `prob` is Hopf-unstable.
Hopf-unstable solutions are real and have exactly two Jacobian eigenvalues with positive real parts,
which are complex conjugates of each other.

`im_tol` : an absolute threshold to distinguish real/complex numbers.
"""
function is_Hopf_unstable(soln::StateDict, res::Result)
    return _is_Hopf_unstable(collect(values(soln)), res.jacobian)
end

"""
$(TYPEDSIGNATURES)

Returns a function that checks if a solution is Hopf-unstable.
"""
function _is_Hopf_unstable(res::Result)
    return _isit(soln) = _is_Hopf_unstable(soln, res.jacobian)
end

"""
$(TYPEDSIGNATURES)

Returns true if the solution `soln` is Hopf-unstable.
"""
function _is_Hopf_unstable(soln, J)
    _is_physical(soln) || return false  # the solution is unphysical anyway
    λs = eigvals(J(soln))
    unstable = filter(x -> real(x) > 0, λs)
    (
        length(unstable) == 2 &&
        abs(conj(unstable[1]) - unstable[2]) < IM_TOL &&
        return true
    ) || return false
    return all(x -> real(x) < 0, λs)
end

"""
$(TYPEDSIGNATURES)

Create binary classification of the solutions, where each solution point receives an identifier based
on its permutation of stable branches (allows to distinguish between different phases, which may
have the same number of stable solutions). It works by converting each bitstring
`[is_stable(solution_1), is_stable(solution_2), ...,]` into unique labels.
"""
function classify_binaries!(res::Result)
    # mapping of binary string (stable or not) for each solution set to integer.
    # Sensitive to ordering!
    bin_label = bitarr_to_int.(clean_bitstrings(res))

    # renormalize labels with numbers from 1 to length(unique(label))
    for (idx, el) in enumerate(unique(bin_label))
        bin_label[findall(x -> x == el, bin_label)] .= idx
    end
    return res.binary_labels .= bin_label
end

"""
$(TYPEDSIGNATURES)

Removes unphysical solutions from the bitstrings in `res`.
"""
function clean_bitstrings(res::Result)
    return [
        [el for el in bit_string[phys_string]] for
        (bit_string, phys_string) in zip(res.classes["stable"], res.classes["physical"])
    ] #remove unphysical solutions from strings
end; #remove unphysical solutions from strings

"""
$(TYPEDSIGNATURES)

Converts a bit array to an integer.
"""
function bitarr_to_int(arr)
    return sum(arr .* (2 .^ collect((length(arr) - 1):-1:0)))
end

"""
$(TYPEDSIGNATURES)

Removes all solution branches from `res` where NONE of the solution falls into `class`.
Typically used to filter out unphysical solutions to prevent huge file sizes.
"""
function filter_result!(res::Result, class::String)
    bools = [any(getindex.(res.classes[class], i)) for i in 1:branch_count(res)]
    res.solutions = [s[bools] for s in res.solutions]
    for c in keys(res.classes)
        res.classes[c] = [s[bools] for s in res.classes[c]]
    end
end

"""
$(TYPEDSIGNATURES)

Classifies solutions in `result` using default classes: "physical", "stable, "Hopf".
"""
function _classify_default!(result)
    classify_solutions!(result, _is_physical, "physical")
    classify_solutions!(result, _is_stable(result), "stable")
    classify_solutions!(result, _is_Hopf_unstable(result), "Hopf")
    order_branches!(result, ["physical", "stable"]) # shuffle the branches to have relevant ones first
    classify_binaries!(result) # assign binaries to solutions depending on which branches are stable
    return nothing
end
