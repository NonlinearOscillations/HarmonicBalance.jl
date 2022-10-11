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
function classify_solutions!(res::Result, condition::String, name::String; physical=true)
    values = classify_solutions(res, condition; physical=physical)
    res.classes[name] = values
end


function classify_solutions(res::Result, condition::String; physical=true)
    expr = Num(eval(Meta.parse(condition)))
    function cond_func(s::OrderedDict, res)
        physical && !is_physical(s, res) && return false
        s = [key => real(s[key]) for key in keys(s)] # make values real
        Bool(substitute_all(expr, s).val)
    end
    classify_solutions(res, cond_func)
end


#   Classify solutions where for `f` which is a function accepting a solution dictionary
#   specifies all params and variables).
function classify_solutions!(res::Result, f::Function, name::String)
    values = classify_solutions(res, f)
    res.classes[name] = values
end


"Return an array of booleans classifying the solution in `res`
according to `f` (`f` takes a solution dictionary, return a boolean)"
function classify_solutions(res::Result, f::Function)
    values = similar(res.solutions, BitVector)
    for (idx, soln) in enumerate(res.solutions)
        values[idx] = [f(get_single_solution(res, index=idx, branch=b), res) for b in 1:length(soln)]
    end
    values
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
function is_physical(soln::StateDict, res::Result; im_tol=im_tol)
    var_values = [getindex(soln, v) for v in res.problem.variables]
    return all(abs.(imag.(var_values)).<im_tol) && any(isnan.(var_values).==false)
end


"""
$(TYPEDSIGNATURES)
Returns true if the solution `soln` of the Result `res` is stable.
Stable solutions are real and have all Jacobian eigenvalues Re[λ] <= 0.
`im_tol` : an absolute threshold to distinguish real/complex numbers.
`rel_tol`: Re(λ) considered <=0 if real.(λ) < rel_tol*abs(λmax)
"""
function is_stable(soln::StateDict, res::Result; im_tol=im_tol, rel_tol=1E-10)
    is_physical(soln, res ,im_tol=im_tol) || return false  # the solution is unphysical anyway
    λs = eigvals(real.(res.jacobian(soln)))
    return all([real.(λs) .< rel_tol*maximum(abs.(λs))]...)
end


"""
$(TYPEDSIGNATURES)
Returns true if the solution `soln` of the problem `prob` is Hopf-unstable.
Hopf-unstable solutions are real and have exactly two Jacobian eigenvalues with +ve real parts, which
are complex conjugates of each other.
`im_tol` : an absolute threshold to distinguish real/complex numbers.
"""
function is_Hopf_unstable(soln::StateDict, res::Result; im_tol=im_tol)
    is_physical(soln, res, im_tol=im_tol) || return false  # the solution is unphysical anyway
    J = res.jacobian(soln)
    unstable = filter(x -> real(x) > 0, eigvals(J))
    (length(unstable) == 2 && abs(conj(unstable[1]) - unstable[2]) < im_tol && return true) || return false
    return all([ real.(eigvals(J)) .< 0]...)
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

