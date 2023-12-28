export get_steady_states
export get_single_solution


"""
$(TYPEDSIGNATURES)
Return an ordered dictionary specifying all variables and parameters of the solution
in `result` on `branch` at the position `index`.
"""
function get_single_solution(res::Result; branch::Int64, index)::OrderedDict{Num, ComplexF64}

    # check if the dimensionality of index matches the solutions
    if length(size(res.solutions)) !== length(index)
        # if index is a number, use linear indexing
        index = length(index) == 1 ? CartesianIndices(res.solutions)[index] : error("Index ", index, " undefined for a solution of size ", size(res.solutions))
    else
        index = CartesianIndex(index)
    end

    vars = OrderedDict(zip(res.problem.variables, res.solutions[index][branch]))

    # collect the swept parameters required for this call
    swept_params = OrderedDict( key => res.swept_parameters[key][index[i]] for (i,key) in enumerate(keys(res.swept_parameters)))
    full_solution = merge(vars, swept_params, res.fixed_parameters)

    return OrderedDict(zip(keys(full_solution), ComplexF64.(values(full_solution))))
end


get_single_solution(res::Result, index) = [get_single_solution(res, index=index, branch = b) for b in 1:length(res.solutions[1])]


"""
    get_steady_states(prob::Problem,
                        swept_parameters::ParameterRange,
                        fixed_parameters::ParameterList;
                        random_warmup=true,
                        threading = Threads.nthreads() > 1,
                        show_progress=true,
                        sorting="nearest")

Solves `prob` over the ranges specified by `swept_parameters`, keeping `fixed_parameters` constant.
`swept_parameters` accepts pairs mapping symbolic variables to arrays or `LinRange`.
`fixed_parameters` accepts pairs mapping symbolic variables to numbers.

Keyword arguments
- `random_warmup`: If `true`, a problem similar to `prob` but with random complex parameters is first solved to find all non-singular paths. The subsequent tracking to find results for all swept_parameters is then much faster than the initial solving. If `random_warmup=false`, each parameter point is solved separately by tracking the maximum number of paths (employs a total degree homotopy).
This takes far longer but can be more reliable.
- `threading`: If `true`, multithreaded support is activated. The number of available threads is set by the environment variable `JULIA_NUM_THREADS`.
- `sorting`: the method used by `sort_solutions` to get continuous solutions branches.  The current options are `"hilbert"` (1D sorting along a Hilbert curve), `"nearest"` (nearest-neighbor sorting) and `"none"`.
- `show_progress`: Indicate whether a progress bar should be displayed.

Example: solving a simple harmonic oscillator ``m \\ddot{x} + γ \\dot{x} + ω_0^2 x = F \\cos(ωt)``
to obtain the response as a function of ``ω``
```julia-repl
# having obtained a Problem object, let's find steady states
julia> range = ParameterRange(ω => LinRange(0.8,1.2,100) ) # 100 parameter sets to solve
julia> fixed = ParameterList(m => 1, γ => 0.01, F => 0.5, ω_0 => 1)
julia> get_steady_states(problem, range, fixed)

A steady state result for 100 parameter points

    Solution branches:   1
       of which real:    1
       of which stable:  1

    Classes: stable, physical, Hopf, binary_labels

```

It is also possible to create multi-dimensional solutions plots.
```julia-repl
# The swept parameters take precedence over fixed -> use the same fixed
julia> range = ParameterRange(ω => LinRange(0.8,1.2,100), F => LinRange(0.1,1.0,10) ) # 100x10 parameter sets

# The swept parameters take precedence over fixed -> the F in fixed is now ignored
julia> get_steady_states(problem, range, fixed)

A steady state result for 1000 parameter points

    Solution branches:   1
       of which real:    1
       of which stable:  1

    Classes: stable, physical, Hopf, binary_labels
```

"""
function get_steady_states(prob::Problem, swept_parameters::ParameterRange, fixed_parameters::ParameterList; random_warmup=true, threading = Threads.nthreads() > 1, show_progress=true, sorting="nearest", classify_default=true)
    # make sure the variables are in our namespace to make them accessible later
    declare_variable.(string.(cat(prob.parameters, prob.variables, dims=1)))

    # prepare an array of vectors, each representing one set of input parameters
    # an n-dimensional sweep uses an n-dimensional array
    unique_fixed = filter_duplicate_parameters(swept_parameters, fixed_parameters)

    any([in(var_name(var), var_name.([keys(fixed_parameters)..., keys(swept_parameters)...])) for var in get_variables(prob)]) && error("Cannot fix one of the variables!")

    input_array = _prepare_input_params(prob, swept_parameters, unique_fixed)
    # feed the array into HomotopyContinuation, get back an similar array of solutions
    raw = _get_raw_solution(prob, input_array, sweep=swept_parameters, random_warmup=random_warmup, threading=threading, show_progress=show_progress)

    # extract all the information we need from results
    #rounded_solutions = unique_points.(HomotopyContinuation.solutions.(getindex.(raw, 1)); metric = EuclideanNorm(), atol=1E-14, rtol=1E-8)
    rounded_solutions = HomotopyContinuation.solutions.(getindex.(raw, 1));
    all(isempty.(rounded_solutions)) ? error("No solutions found!") : nothing
    solutions = pad_solutions(rounded_solutions)

    compiled_J = _compile_Jacobian(prob, swept_parameters, unique_fixed)

    result = Result(solutions, swept_parameters, unique_fixed, prob, Dict(), compiled_J) # a "raw" solution struct

    sort_solutions!(result, sorting=sorting, show_progress=show_progress) # sort into branches
    classify_default ? _classify_default!(result) : nothing

    return result

end


function _classify_default!(result)
    classify_solutions!(result, _is_physical, "physical")
    classify_solutions!(result, _is_stable(result), "stable")
    classify_solutions!(result, _is_Hopf_unstable(result), "Hopf")
    order_branches!(result, ["physical", "stable"]) # shuffle the branches to have relevant ones first
    classify_binaries!(result) # assign binaries to solutions depending on which branches are stable
end


get_steady_states(p::Problem, swept, fixed; kwargs...) = get_steady_states(p, ParameterRange(swept), ParameterList(fixed); kwargs...)
get_steady_states(eom::HarmonicEquation, swept, fixed; kwargs...) = get_steady_states(Problem(eom), swept, fixed; kwargs...)
get_steady_states(p, pairs; kwargs...) = get_steady_states(p, filter( x->length(x[2]) > 1, pairs), filter(x -> length(x[2]) == 1, pairs); kwargs...)


""" Compile the Jacobian from `prob`, inserting `fixed_parameters`.
    Returns a function that takes a dictionary of variables and `swept_parameters` to give the Jacobian."""
function _compile_Jacobian(prob::Problem, swept_parameters::ParameterRange, fixed_parameters::ParameterList)
    J_variables = cat(prob.variables, collect(keys(swept_parameters)), dims=1)

    if prob.jacobian isa Matrix
        compiled_J = compile_matrix(prob.jacobian, J_variables, fixed_parameters)
    else
        compiled_J = prob.jacobian # leave implicit Jacobian as is
    end
    compiled_J
end

"""
Take a matrix containing symbolic variables `variables` and keys of `fixed_parameters`.
Substitute the values according to `fixed_parameters` and compile into a function that takes numerical arguments
    in the order set in `variables`.
"""
function compile_matrix(matrix, variables, fixed_parameters)
    J = substitute_all(matrix, fixed_parameters)
    matrix = build_function(J, variables)
    matrix = eval(matrix[1]) # compiled allocating function, see Symbolics manual
    m(vals::Vector) = matrix(vals)
    m(s::OrderedDict) = m([s[var] for var in variables]) # for the UI
    return m
end


"Find a branch order according `classification`. Place branches where true occurs earlier first."
function find_branch_order(classification::Vector{BitVector})
    branches = [getindex.(classification, k) for k in 1:length(classification[1])] # array of branches
    indices = replace(findfirst.(branches), nothing => Inf)
    negative = findall(x -> x == Inf, indices) # branches not true anywhere - leave out
    order = setdiff(sortperm(indices), negative)
end

find_branch_order(classification::Array) = collect(1:length(classification[1])) # no ordering for >1D

"Order the solution branches in `res` such that close classified positively by `classes` are first."
function order_branches!(res::Result, classes::Vector{String})
    for class in classes
        order_branches!(res, find_branch_order(res.classes[class]))
    end
end

order_branches!(res::Result, class::String) = order_branches!(res, [class])

"Reorder the solutions in `res` to match the index permutation `order`."
function order_branches!(res::Result, order::Vector{Int64})
    res.solutions = _reorder_nested(res.solutions, order)
    for key in keys(res.classes)
        res.classes[key] = _reorder_nested(res.classes[key], order)
    end
end

"Reorder EACH ELEMENT of `a` to match the index permutation `order`. If length(order) < length(array), the remanining positions are kept."
function _reorder_nested(a::Array, order::Vector{Int64})
    a[1] isa Union{Array, BitVector} || return a
    order = length(order) == length(a) ? order : vcat(order, setdiff(1:length(a[1]), order)) # pad if needed
    new_array = [el[order] for el in a]
end


"prepares an input vector to be parsed to the 2D phase diagram with parameters to sweep and kwargs"
function _prepare_input_params(prob::Problem, sweeps::ParameterRange, fixed_parameters::ParameterList)
    # sweeping takes precedence over fixed_parameters
    fixed_parameters = filter_duplicate_parameters(sweeps, fixed_parameters)

    # fixed order of parameters
    all_keys = cat(collect(keys(sweeps)), collect(keys(fixed_parameters)), dims=1)
    # the order of parameters we have now does not correspond to that in prob!
    # get the order from prob and construct a permutation to rearrange our parameters
    error = ArgumentError("Some input parameters are missing or appear more than once!")
    permutation = try
         p = [findall(x->isequal(x, par), all_keys) for par in prob.parameters] # find the matching position of each parameter
         all((length.(p)) .== 1) || throw(error) # some parameter exists more than twice!
         p = getindex.(p, 1) # all exist once -> flatten
         isequal(all_keys[getindex.(p, 1)], prob.parameters) || throw(error) # parameters sorted wrong!
         #isempty(setdiff(all_keys, prob.parameters)) || throw(error) # extra parameters present!
         p
    catch
        throw(error)
    end

    param_ranges = collect(values(sweeps)) # Vector of the sweep LinRanges
    input_array = collect(Iterators.product(param_ranges..., values(fixed_parameters)...)) # array of all permutations (fixed_params do not change)

    # order each parameter vector to match the order in prob
    input_array = getindex.(input_array, [permutation])
    # HC wants arrays, not tuples
    return tuple_to_vector.(input_array)
end


"Remove occurrences of `sweeps` elements from `fixed_parameters`."
function filter_duplicate_parameters(sweeps, fixed_parameters)
    new_params = copy(fixed_parameters)
    for par in keys(sweeps)
        delete!(new_params, par)
    end
    return new_params
end


"Uses HomotopyContinuation to solve `problem` at specified `parameter_values`."
function _get_raw_solution(problem::Problem, parameter_values; sweep=[], random_warmup=false, threading=false, show_progress=true)
    # HomotopyContinuation accepts 1D arrays of parameter sets
    params_1D = reshape(parameter_values, :, 1)

    if random_warmup && !isempty(sweep)
        complex_pert = [1E-2 * issubset(p, keys(sweep))*randn(ComplexF64) for p in problem.parameters] # complex perturbation of the warmup parameters
        warmup_parameters = params_1D[Int(round(length(params_1D)/2))] .* (ones(length(params_1D[1])) + complex_pert)
        warmup_solution = HomotopyContinuation.solve(problem.system,  target_parameters=warmup_parameters, threading=threading, show_progress=show_progress)
        result_full = HomotopyContinuation.solve(problem.system, HomotopyContinuation.solutions(warmup_solution), start_parameters=warmup_parameters, target_parameters=parameter_values, threading=threading, show_progress=show_progress)
    else
        result_full = Array{Vector{Any}, 1}(undef, length(parameter_values))
        if show_progress
            bar = Progress(length(parameter_values), 1, "Solving via total degree homotopy ...", 50)
        end
        for i in 1:length(parameter_values) # do NOT thread this
            p = parameter_values[i]
            show_progress ? next!(bar) : nothing
            result_full[i] = [HomotopyContinuation.solve(problem.system, start_system=:total_degree, target_parameters=p, threading=threading, show_progress=false), p]
        end
    end

    # reshape back to the original shape
    return reshape(result_full, size(parameter_values)...)
end


"Add `padding_value` to `solutions` in case their number changes in parameter space."
function pad_solutions(solutions::Array{Vector{Vector{ComplexF64}}}; padding_value=NaN)
    Ls    = length.(solutions)
    nvars = length(solutions[1][1]) # number of variables
    max_N = maximum(Ls) # length to be fixed
    padded_solutions = deepcopy(solutions)
    for (i,s) in enumerate(solutions)
        if Ls[i]<max_N
            padded_solutions[i] = vcat(solutions[i],[padding_value*ones(nvars) for k in 1:(max_N-length(s))])
        end
    end
    return padded_solutions
end


tuple_to_vector(t::Tuple) = [i for i in t]


function newton(prob::Problem, soln::OrderedDict)

    vars = _convert_or_zero.(substitute_all(prob.variables, soln))
    pars = _convert_or_zero.(substitute_all(prob.parameters, soln))

    HomotopyContinuation.newton(prob.system, vars, pars)
end


"""
    newton(res::Result, soln::OrderedDict)
    newton(res::Result; branch, index)

Run a newton solver on `prob::Problem` starting from the solution `soln` (indexable by `branch` and `index`).
Any variables/parameters not present in `soln` are set to zero.
"""
newton(res::Result, soln::OrderedDict) = newton(res.problem, soln)
newton(res::Result; branch, index) = newton(res, res[index][branch])


function _convert_or_zero(x, t=ComplexF64)
    try
        convert(t, x)
    catch ArgumentError
        @warn string(x) * " not supplied: setting to zero"
        return 0
    end
end
