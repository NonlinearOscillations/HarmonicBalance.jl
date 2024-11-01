"""
    get_steady_states(prob::Problem,
                        swept_parameters::ParameterRange,
                        fixed_parameters::ParameterList;
                        method=:warmup,
                        threading = Threads.nthreads() > 1,
                        show_progress=true,
                        sorting="nearest")

Solves `prob` over the ranges specified by `swept_parameters`, keeping `fixed_parameters` constant.
`swept_parameters` accepts pairs mapping symbolic variables to arrays or `LinRange`.
`fixed_parameters` accepts pairs mapping symbolic variables to numbers.

Keyword arguments
- `method`: If `:warmup` (default), a problem similar to `prob` but with random complex parameters is first solved to find all non-singular paths. The subsequent tracking to find results for all `swept_parameters` is then much faster than the initial solving. If `method=:total_degree`, each parameter point is solved separately by tracking the maximum number of paths (employs a total degree homotopy).
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
function get_steady_states(
    prob::Problem,
    method::HarmonicBalanceMethod,
    swept_parameters::ParameterRange,
    fixed_parameters::ParameterList;
    show_progress=true,
    sorting="nearest",
    classify_default=true,
)
    Random.seed!(seed(method))
    # make sure the variables are in our namespace to make them accessible later
    declare_variable.(string.(cat(prob.parameters, prob.variables; dims=1)))

    variable_names = var_name.([keys(fixed_parameters)..., keys(swept_parameters)...])
    any([var_name(var) ∈ variable_names for var in get_variables(prob)]) && error("Cannot fix one of the variables!")

    unique_fixed, input_array = _prepare_input_params(
        prob, swept_parameters, fixed_parameters
    )
    solutions = get_solutions(prob, method, input_array; show_progress=show_progress)

    compiled_J = _compile_Jacobian(prob, swept_parameters, unique_fixed)

    result = Result(
        solutions, swept_parameters, unique_fixed, prob, Dict(), compiled_J, seed(method)
    )

    if sorting != "no_sorting"
        sort_solutions!(result; sorting=sorting, show_progress=show_progress)
    end
    classify_default ? _classify_default!(result) : nothing

    return result
end

function get_steady_states(
    p::Problem, method::HarmonicBalanceMethod, swept, fixed; kwargs...
)
    return get_steady_states(
        p, method, ParameterRange(swept), ParameterList(fixed); kwargs...
    )
end
function get_steady_states(
    eom::HarmonicEquation, method::HarmonicBalanceMethod, swept, fixed; kwargs...
)
    return get_steady_states(Problem(eom), method, swept, fixed; kwargs...)
end
function get_steady_states(eom::HarmonicEquation, pairs::Dict; kwargs...)
    swept = filter(x -> length(x[2]) > 1, pairs)
    fixed = filter(x -> length(x[2]) == 1, pairs)
    return get_steady_states(eom, swept, fixed; kwargs...)
end
function get_steady_states(
    eom::HarmonicEquation, method::HarmonicBalanceMethod, pairs::Dict; kwargs...
)
    swept = filter(x -> length(x[2]) > 1, pairs)
    fixed = filter(x -> length(x[2]) == 1, pairs)
    return get_steady_states(eom, method, swept, fixed; kwargs...)
end
function get_steady_states(eom::HarmonicEquation, swept, fixed; kwargs...)
    return get_steady_states(Problem(eom), WarmUp(), swept, fixed; kwargs...)
end

function get_solutions(prob, method, input_array; show_progress)
    raw = _get_raw_solution(prob, method, input_array; show_progress=show_progress)

    solutions = HC.solutions.(getindex.(raw, 1))
    if all(isempty.(solutions))
        @warn "No solutions found!"
        return solutions
    else
        pad_solutions(solutions)
    end
end

"""
Take a matrix containing symbolic variables `variables` and keys of `fixed_parameters`.
Substitute the values according to `fixed_parameters` and compile into a function that takes numerical arguments
    in the order set in `variables`.
"""
function compile_matrix(mat, variables; rules=Dict(), postproc=x -> x)
    J = substitute_all.(mat, Ref(rules))
    matrix = Symbolics.build_function(J, variables)
    matrix = eval(matrix[1]) # compiled allocating function, see Symbolics manual
    m(vals::Vector) = postproc(matrix(vals))
    m(s::OrderedDict) = m([s[var] for var in variables]) # for the UI
    return m
end

"Reorder EACH ELEMENT of `a` to match the index permutation `order`. If length(order) < length(array), the remanining positions are kept."
function _reorder_nested(a::Array, order::Vector{Int64})
    a[1] isa Union{Array,BitVector} || return a
    order = length(order) == length(a) ? order : vcat(order, setdiff(1:length(a[1]), order)) # pad if needed
    return new_array = [el[order] for el in a]
end

"prepares an input vector to be parsed to the 2D phase diagram with parameters to sweep and kwargs"
function _prepare_input_params(
    prob::Problem, sweeps::ParameterRange, fixed_parameters::ParameterList
)
    # sweeping takes precedence over fixed_parameters
    unique_fixed = filter_duplicate_parameters(sweeps, fixed_parameters)

    # fixed order of parameters
    all_keys = cat(collect(keys(sweeps)), collect(keys(fixed_parameters)); dims=1)
    # the order of parameters we have now does not correspond to that in prob!
    # get the order from prob and construct a permutation to rearrange our parameters
    error = ArgumentError("Some input parameters are missing or appear more than once!")
    permutation = try
        p = [findall(x -> isequal(x, par), all_keys) for par in prob.parameters] # find the matching position of each parameter
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
    return unique_fixed, tuple_to_vector.(input_array)
end

"Remove occurrences of `sweeps` elements from `fixed_parameters`."
function filter_duplicate_parameters(sweeps, fixed_parameters)
    new_params = copy(fixed_parameters)
    for par in keys(sweeps)
        delete!(new_params, par)
    end
    return new_params
end

"A random warmup solution is computed to use as `start_parameters` in the homotopy."
function _solve_warmup(problem::Problem, method::WarmUp, params; show_progress)
    # complex perturbation of the warmup parameters
    options = alg_specific_options(method)
    l = length(params[1])

    perturbation = ones(l) + options[:perturbation_size] * randn(ComplexF64, l)
    warmup_parameters = params[options[:index]] .* perturbation

    warmup_solution = HC.solve(
        problem.system;
        start_system=:total_degree,
        target_parameters=warmup_parameters,
        show_progress=show_progress,
        alg_default_options(method)...,
    )
    return warmup_parameters, warmup_solution
end

"Uses HomotopyContinuation to solve `problem` at specified `parameter_values`."
function _get_raw_solution(
    problem::Problem, method::WarmUp, parameter_values; show_progress
)
    warmup_parameters, warmup_solution = _solve_warmup(
        problem, method, parameter_values; show_progress=show_progress
    )
    result_full = HC.solve(
        problem.system,
        HC.solutions(warmup_solution);
        start_parameters=warmup_parameters,
        target_parameters=parameter_values,
        show_progress=show_progress,
        alg_default_options(method)...,
    )

    return reshape(result_full, size(parameter_values)...)
end

"Uses HomotopyContinuation to solve `problem` at specified `parameter_values`."
function _get_raw_solution(
    problem::Problem, method::Union{TotalDegree,Polyhedral}, parameter_values; show_progress
)
    result_full = Array{Vector{Any},1}(undef, length(parameter_values))
    if show_progress
        bar = Progress(
            length(parameter_values);
            dt=1,
            desc="Solving via $method homotopy ...",
            barlen=50,
        )
    end
    for i in eachindex(parameter_values) # do NOT thread this
        p = parameter_values[i]
        show_progress ? ProgressMeter.next!(bar) : nothing
        result_full[i] = [
            HC.solve(
                problem.system;
                start_system=method_symbol(method),
                target_parameters=p,
                show_progress=false,
                alg_default_options(method)...,
                alg_specific_options(method)...,
            ),
            p,
        ]
    end

    # reshape back to the original shape
    return reshape(result_full, size(parameter_values)...)
end

"Add `padding_value` to `solutions` in case their number changes in parameter space."
function pad_solutions(solutions::Array{Vector{Vector{ComplexF64}}}; padding_value=NaN)
    Ls = length.(solutions)
    nvars = length(solutions[1][1]) # number of variables
    max_N = maximum(Ls) # length to be fixed
    padded_solutions = deepcopy(solutions)
    for (i, s) in enumerate(solutions)
        if Ls[i] < max_N
            padded_solutions[i] = vcat(
                solutions[i], [padding_value * ones(nvars) for k in 1:(max_N - length(s))]
            )
        end
    end
    return padded_solutions
end

function newton(prob::Problem, soln::OrderedDict)
    vars = _convert_or_zero.(substitute_all(prob.variables, soln))
    pars = _convert_or_zero.(substitute_all(prob.parameters, soln))

    return HC.newton(prob.system, vars, pars)
end

"""
    newton(res::Result, soln::OrderedDict)
    newton(res::Result; branch, index)

Run a newton solver on `prob::Problem` starting from the solution `soln` (indexable by `branch` and `index`).
Any variables/parameters not present in `soln` are set to zero.
"""
newton(res::Result, soln::OrderedDict) = newton(res.problem, soln)
newton(res::Result; branch, index) = newton(res, res[index][branch])
