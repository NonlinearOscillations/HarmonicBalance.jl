"""
    get_steady_states(problem::HarmonicEquation,
                        method::HarmonicBalanceMethod,
                        swept_parameters,
                        fixed_parameters;
                        show_progress=true,
                        sorting="nearest",
                        classify_default=true)

Solves `problem` with the `method` over the ranges specified by `swept_parameters`,
keeping `fixed_parameters` constant.
`swept_parameters` accepts pairs mapping symbolic variables to arrays or ranges.
`fixed_parameters` accepts pairs mapping symbolic variables to numbers.

### Keyword arguments
- `show_progress`: Indicate whether a progress bar should be displayed.
- `sorting`: the method used by `sort_solutions` to get continuous solutions branches.
    The current options are `"hilbert"` (1D sorting along a Hilbert curve), `"nearest"`
    (nearest-neighbor sorting) and `"none"`.
- `classify_default`: If `true`, the solutions will be classified using the default
    classification method.

### Example
solving a simple harmonic oscillator
``m \\ddot{x} + γ \\dot{x} + ω_0^2 x = F \\cos(ωt)`` to obtain the response
as a function of ``ω``
```julia-repl
# having obtained a Problem object, let's find steady states
julia> range = (ω => range(0.8, 1.2, 100) ) # 100 parameter sets to solve
julia> fixed = ParameterList(m => 1, γ => 0.01, F => 0.5, ω_0 => 1)
julia> get_steady_states(problem, range, fixed)

A steady state result for 100 parameter points

    Solution branches:   1
       of which real:    1
       of which stable:  1

    Classes: stable, physical, Hopf, binary_labels

```

It is also possible to perform 2-dimensional sweeps.
```julia-repl
# The swept parameters take precedence over fixed -> use the same fixed
julia> range = (ω => range(0.8,1.2,100), F => range(0.1,1.0,10) )

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
    method::HarmonicBalanceMethod;
    show_progress=true,
    sorting="nearest",
    classify_default=true,
    verbose=false,
)
    Random.seed!(seed(method))
    # make sure the variables are in our namespace to make them accessible later
    declare_variables(prob)

    swept_parameters = prob.swept_parameters
    fixed_parameters = prob.fixed_parameters

    unique_fixed, input_array = _prepare_input_params(
        prob, swept_parameters, fixed_parameters
    )
    verbose && @info "Find solutions"
    solutions = get_solutions(prob, method, input_array; show_progress)

    result = Result(
        solutions,
        swept_parameters,
        unique_fixed,
        prob,
        Dict(),
        zeros(Int64, size(solutions)...),
        prob.jacobian,
        seed(method),
    )

    verbose && @info "Sort solutions"
    if sorting != "no_sorting"
        sort_solutions!(result; sorting, show_progress)
    end

    verbose && @info "Classify solutions"
    classify_default ? _classify_default!(result) : nothing

    return result
end

function get_steady_states(
    eom::HarmonicEquation, method::HarmonicBalanceMethod, swept, fixed; kwargs...
)
    return get_steady_states(
        Problem(eom, OrderedDict(swept), OrderedDict(fixed)), method; kwargs...
    )
end
function get_steady_states(eom::HarmonicEquation, pairs::Union{Dict,OrderedDict}; kwargs...)
    swept = filter(x -> length(x[2]) > 1, pairs)
    fixed = filter(x -> length(x[2]) == 1, pairs)
    return get_steady_states(eom, swept, fixed; kwargs...)
end
function get_steady_states(
    eom::HarmonicEquation,
    method::HarmonicBalanceMethod,
    pairs::Union{Dict,OrderedDict};
    kwargs...,
)
    swept = filter(x -> length(x[2]) > 1, pairs)
    fixed = filter(x -> length(x[2]) == 1, pairs)
    return get_steady_states(eom, method, swept, fixed; kwargs...)
end
function get_steady_states(eom::HarmonicEquation, swept, fixed; kwargs...)
    return get_steady_states(
        Problem(eom, OrderedDict(swept), OrderedDict(fixed)), WarmUp(); kwargs...
    )
end

function get_solutions(prob, method, input_array; show_progress)
    raw = _get_raw_solution(prob, method, input_array; show_progress)

    if all(isempty.(raw))
        @warn "No solutions found!"
        return raw
    else
        pad_solutions(raw)
    end
end

"""
Reorder EACH ELEMENT of `a` to match the index permutation `order`.
If length(order) < length(array), the remanining positions are kept.
"""
function _reorder_nested(a::Array, order::Vector{Int})
    a[1] isa Union{Array,BitVector} || return a
    order = length(order) == length(a) ? order : vcat(order, setdiff(1:length(a[1]), order)) # pad if needed
    return new_array = [el[order] for el in a]
end

"""
prepares an input vector to be parsed to the 2D phase diagram with parameters
to sweep and kwargs
"""
function _prepare_input_params(
    prob::Problem, sweeps::OrderedDict, fixed_parameters::OrderedDict
)
    unique_fixed, permutation = unique_fixed_and_permutations(
        prob, sweeps, fixed_parameters
    )

    input_array = type_stable_parameters(sweeps, fixed_parameters)
    # order each parameter vector to match the order in prob
    input_array = getindex.(input_array, [permutation])
    # HC wants arrays, not tuples
    return unique_fixed, tuple_to_vector.(input_array)
end

"Uses HomotopyContinuation to solve `problem` at specified `parameter_values`."
function _get_raw_solution(
    problem::Problem, method::WarmUp, parameter_values; show_progress
)
    warm_up_method = method.warm_up_method

    example_p = parameter_values[1]
    if isempty(method.start_parameters)
        start_parameters = randn(complex(eltype(example_p)), length(example_p))
    else
        start_parameters = method.start_parameters
    end

    warmup_solution = HC.solve(
        problem.system;
        start_system=method_symbol(warm_up_method),
        target_parameters=start_parameters,
        show_progress,
        alg_default_options(warm_up_method)...,
        alg_specific_options(warm_up_method)...,
    )

    result_full = HC.solve(
        problem.system,
        HC.solutions(warmup_solution);
        start_parameters=start_parameters,
        target_parameters=parameter_values,
        show_progress,
        alg_default_options(method)...,
    )

    raw_solutions = reshape(result_full, size(parameter_values)...)
    raw_solutions = HC.solutions.(getindex.(raw_solutions, 1))
    # cache = HC.NewtonCache(F)
    if method.check_zero
        Threads.@threads for i in eachindex(parameter_values)
            any(is_zero.(raw_solutions[i])) && continue
            p = parameter_values[i]
            zero_root = HC.newton(problem.system, zeros(length(problem.variables)), p)
            if is_zero(zero_root.x)
                push!(raw_solutions[i], zero_root.x)
            end
        end
    end

    return raw_solutions
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
            desc="Solving via $(nameof(typeof(method))) homotopy ...",
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

    raw_solutions = reshape(result_full, size(parameter_values)...)
    raw_solutions = HC.solutions.(getindex.(raw_solutions, 1))

    # reshape back to the original shape
    return raw_solutions
end

"Add `padding_value` to `solutions` in case their number changes in parameter space."
function pad_solutions(solutions::Solutions(T); padding_value=NaN) where {T}
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
    vars = substitute_all(prob.variables, soln)
    pars = substitute_all(prob.parameters, soln)

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
