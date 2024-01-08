export classify_branch
export classify_solutions!

"""
$(TYPEDSIGNATURES)

Creates a solution class in `res` using the inequality `condition` (parsed into Symbolics.jl input).

The new class is labelled with `name` and stored under `res.classes[name]`.

By default, only physical (=real) solutions are classified, `false` is returned for the rest.

# Example
```julia
# solve a previously-defined problem
res = get_steady_states(problem, swept_parameters, fixed_parameters)

# classify, store in result.classes["large_amplitude"]
classify_solutions!(res, "sqrt(u1^2 + v1^2) > 1.0" , "large_amplitude")
```

"""
function classify_solutions!(res::Result, func::Union{String, Function}, name::String; physical=true, kwargs...)
    values = classify_solutions(res, func; physical=physical, kwargs...)
    res.classes[name] = values
end


function classify_solutions(res::Result, func; physical=true, kwargs...)
    func = isa(func, Function) ? func : _build_substituted(func, res)
    if physical
        f_comp(soln; kwargs...) = _is_physical(soln) && func(real.(soln); kwargs...)
        transform_solutions(res, f_comp; kwargs...)

    else
        transform_solutions(res, func; kwargs...)
    end
end


"""
$(TYPEDSIGNATURES)
Returns an array of booleans classifying `branch` in the solutions in `res`
according to `class`.
"""
function classify_branch(res::Result, branch::Int64, class::String)
    branch_values = getindex.(res.classes[class], branch)
end

classify_branch(soln::Result, class::String) = [classify_branch(soln, b, class) for b in 1:length(first(soln.solutions))]


"""
$(TYPEDSIGNATURES)
Returns true if the solution `soln` of the Result `res` is physical (= real number).
`im_tol` : an absolute threshold to distinguish real/complex numbers.
"""
function is_physical(soln::StateDict, res::Result)
    var_values = [getindex(soln, v) for v in res.problem.variables]
    return _is_physical(var_values)
end

_is_physical(soln; im_tol=IM_TOL) = all( x -> !isnan(x) && abs(imag(x)) < im_tol, soln)
_is_physical(res::Result) = classify_solutions(res, _is_physical)


"""
$(TYPEDSIGNATURES)
Returns true if the solution `soln` of the Result `res` is stable.
Stable solutions are real and have all Jacobian eigenvalues Re[λ] <= 0.
`im_tol` : an absolute threshold to distinguish real/complex numbers.
`rel_tol`: Re(λ) considered <=0 if real.(λ) < rel_tol*abs(λmax)
"""
function is_stable(soln::StateDict, res::Result; kwargs...)
    _is_stable(values(soln) |> collect, res.jacobian; kwargs...)
end

function _is_stable(res::Result; kwargs...)
    _isit(soln) = _is_stable(soln, res.jacobian; kwargs...)
end

function _is_stable(soln, J; rel_tol=1E-10)
    _is_physical(soln) || return false
    λs = eigvals(real.(J(soln)))
    scale = maximum(Iterators.map(abs, λs))
    all(x -> real(x) < rel_tol*scale, λs)
end

"""
$(TYPEDSIGNATURES)
Returns true if the solution `soln` of the problem `prob` is Hopf-unstable.
Hopf-unstable solutions are real and have exactly two Jacobian eigenvalues with +ve real parts, which
are complex conjugates of each other.
`im_tol` : an absolute threshold to distinguish real/complex numbers.
"""
function is_Hopf_unstable(soln::StateDict, res::Result)
    _is_Hopf_unstable(values(soln) |> collect, res.jacobian)
end

function _is_Hopf_unstable(res::Result)
    _isit(soln) = _is_Hopf_unstable(soln, res.jacobian)
end

function _is_Hopf_unstable(soln, J)
    _is_physical(soln) || return false  # the solution is unphysical anyway
    λs = eigvals(J(soln))
    unstable = filter(x -> real(x) > 0, λs)
    (length(unstable) == 2 && abs(conj(unstable[1]) - unstable[2]) < IM_TOL && return true) || return false
    return all(x -> real(x) < 0, λs)
end


"""
$(TYPEDSIGNATURES)
Create binary classification of the solutions, such that each solution points receives an identifier based
on its permutation of stable branches (allows to distinguish between different phases, which may have the same number
of stable solutions). It works by converting each bistring `[is_stable(solution_1),is_stable(solution_2),...,]` into unique labels.
"""
function classify_binaries!(res::Result)
    bin_label  = bitarr_to_int.(clean_bitstrings(res)) #mapping of binary string (stable or not) for each solution set to integer. Sensitive to ordering!
    #renormalize labels with numbers from 1 to length(unique(label))
    for (idx,el) in enumerate(unique(bin_label))
        bin_label[findall(x->x==el, bin_label)] .= idx
    end
    res.classes["binary_labels"] = bin_label
end

clean_bitstrings(res::Result) = [[el for el in bit_string[phys_string]]
for (bit_string,phys_string) in zip(res.classes["stable"],res.classes["physical"])]; #remove unphysical solutions from strings

function bitarr_to_int(arr)
    return sum(arr .* (2 .^ collect(length(arr)-1:-1:0)))
end


"""
$(TYPEDSIGNATURES)
Removes all solution branches from `res` where NONE of the solution falls into `class`.
Typically used to filter out unphysical solutions to prevent huge file sizes.
"""
function filter_result!(res::Result, class::String)
    bools = [any(getindex.(res.classes[class],i)) for i in 1:length(res[1])]
    res.solutions = [s[bools] for s in res.solutions]
    for c in filter(x -> x != "binary_labels", keys(res.classes)) # binary_labels stores one Int per parameter set, ignore here
        res.classes[c] = [s[bools] for s in res.classes[c]]
    end
end
