"""
$(TYPEDSIGNATURES)
Return an ordered dictionary specifying all variables and parameters of the solution
in `result` on `branch` at the position `index`.
"""
function get_single_solution(res::Result; branch::Int64, index)::OrderedDict{Num,ComplexF64}

    # check if the dimensionality of index matches the solutions
    if length(size(res.solutions)) !== length(index)
        # if index is a number, use linear indexing
        index = if length(index) == 1
            CartesianIndices(res.solutions)[index]
        else
            error("Index ", index, " undefined for a solution of size ", size(res.solutions))
        end
    else
        index = CartesianIndex(index)
    end

    vars = OrderedDict(zip(res.problem.variables, res.solutions[index][branch]))

    # collect the swept parameters required for this call
    swept_params = OrderedDict(
        key => res.swept_parameters[key][index[i]] for
        (i, key) in enumerate(keys(res.swept_parameters))
    )
    full_solution = merge(vars, swept_params, res.fixed_parameters)

    return OrderedDict(zip(keys(full_solution), ComplexF64.(values(full_solution))))
end

function get_single_solution(res::Result, index)
    return [
        get_single_solution(res; index=index, branch=b) for b in 1:length(res.solutions[1])
    ]
end

"""
$(TYPEDSIGNATURES)

Takes a `Result` object and a string `f` representing a Symbolics.jl expression.
Returns an array with the values of `f` evaluated for the respective solutions.
Additional substitution rules can be specified in `rules` in the format `("a" => val)` or `(a => val)`
"""
function transform_solutions(res::Result, func; branches=1:branch_count(res), realify=false)
    # preallocate an array for the numerical values, rewrite parts of it
    # when looping through the solutions
    pars = collect(values(res.swept_parameters))
    n_vars = length(get_variables(res))
    n_pars = length(pars)

    vtype = if isa(Base.invokelatest(func, rand(Float64, n_vars + n_pars)), Bool)
        BitVector
    else
        Vector{ComplexF64}
    end
    transformed = _similar(vtype, res; branches=branches)
    f = realify ? v -> real.(v) : identity

    batches = Iterators.partition(
        CartesianIndices(res.solutions),
        ceil(Int, length(res.solutions) / Threads.nthreads()),
    )
    Threads.@threads for batch in collect(batches)
        _vals = Vector{ComplexF64}(undef, n_vars + n_pars)
        for idx in batch
            for i in 1:length(idx) # param values are common to all branches
                _vals[end - n_pars + i] = pars[i][idx[i]]
            end
            for (k, branch) in enumerate(branches)
                _vals[1:n_vars] .= res.solutions[idx][branch]
                transformed[idx][k] = Base.invokelatest(func, f(_vals)) # beware, func may be mutating
            end
        end
    end
    return transformed
end

function transform_solutions(res::Result, f::String; rules=Dict(), kwargs...)
    # a string is used as input
    # a macro would not "see" the user's namespace while the user's namespace does not "see" the variables
    func = _build_substituted(f, res; rules=rules)
    return transform_solutions(res, func; kwargs...)
end

function transform_solutions(res::Result, fs::Vector{String}; kwargs...)
    return [transform_solutions(res, f; kwargs...) for f in fs]
end

# a simplified version meant to work with arrays of solutions
# cannot parse parameter values -- meant for time-dependent results
function transform_solutions(soln::Vector, f::String, harm_eq::HarmonicEquation)
    vars = _remove_brackets(get_variables(harm_eq))
    transformed = Vector{ComplexF64}(undef, length(soln))

    # parse the input with Symbolics
    expr = _parse_expression(f)

    rule(u) = Dict(zip(vars, u))

    transformed = map(x -> Symbolics.unwrap(substitute_all(expr, rule(x))), soln)
    return convert(typeof(soln[1]), transformed)
end

""" Parse `expr` into a Symbolics.jl expression, substitute with `rules` and build a function taking free_symbols """
function _build_substituted(expr::String, rules, free_symbols)
    subbed = substitute_all(_parse_expression(expr), rules)
    comp_func = Symbolics.build_function(subbed, free_symbols)

    return eval(comp_func)
end

""" Parse `expr` into a Symbolics.jl expression, substituting the fixed parameters of `res`
The resulting function takes in the values of the variables and swept parameters. """
function _build_substituted(expr, res::Result; rules=Dict())

    # define variables in rules in this namespace
    new_keys = declare_variable.(string.(keys(Dict(rules))))
    fixed_subs = merge(res.fixed_parameters, Dict(zip(new_keys, values(Dict(rules)))))

    return _build_substituted(expr, fixed_subs, _free_symbols(res))
end

function _similar(type, res::Result; branches=1:branch_count(res))
    return [type(undef, length(branches)) for k in res.solutions]
end

"""
$(TYPEDSIGNATURES)

Return an array of bools to mark solutions in `res` which fall into `classes` but not `not_classes`.
Only `branches` are considered.
"""
function _get_mask(res, classes, not_classes=[]; branches=1:branch_count(res))
    classes == "all" && return fill(trues(length(branches)), size(res.solutions))
    bools = vcat(
        [res.classes[c] for c in _str_to_vec(classes)],
        [map(.!, res.classes[c]) for c in _str_to_vec(not_classes)],
    )
    #m = map( x -> [getindex(x, b) for b in [branches...]], map(.*, bools...))

    return m = map(x -> x[[branches...]], map(.*, bools...))
end

"""
$(TYPEDSIGNATURES)

Go over a solution and an equally-sized array (a "mask") of booleans.
true  -> solution unchanged
false -> changed to NaN (omitted from plotting)
"""
function _apply_mask(solns::Array{Vector{ComplexF64}}, booleans)
    factors = replace.(booleans, 0 => NaN)
    return map(.*, solns, factors)
end
function _apply_mask(solns::Vector{Vector{Vector{ComplexF64}}}, booleans)
    Nan_vector = NaN .* similar(solns[1][1])
    new_solns = [
        [booleans[i][j] ? solns[i][j] : Nan_vector for j in eachindex(solns[i])] for
        i in eachindex(solns)
    ]
    return new_solns
end

###
# TRANSFORMATIONS TO THE LAB frame
###

"""
    to_lab_frame(res::Result, x::Num, times; index::Int, branch::Int)
    to_lab_frame(soln::OrderedDict, res::Result, nat_var::Num, times)

Transform a solution into the lab frame (i.e., invert the harmonic ansatz) for the natural variable `x` for `times`. You can also compute the velocity by passing `d(x,t)` for the natural variable `x`.
Either extract the solution from `res::Result` by `index` and `branch` or input `soln::OrderedDict` explicitly.
"""
function to_lab_frame(soln::OrderedDict, res::Result, nat_var::Num, times)
    count_derivatives(nat_var) > 2 &&
        throw(ArgumentError("Only the first derivative is supported"))

    var = if Symbolics.is_derivative(unwrap(nat_var))
        wrap(first(Symbolics.arguments(unwrap(nat_var))))
    else
        nat_var
    end
    vars = filter(x -> isequal(x.natural_variable, var), res.problem.eom.variables)

    return if Symbolics.is_derivative(unwrap(nat_var))
        _to_lab_frame_velocity(soln, vars, times)
    else
        _to_lab_frame(soln, vars, times)
    end
end
function to_lab_frame(res::Result, nat_var::Num, times; index::Int, branch::Int)
    return to_lab_frame(res[index][branch], res, nat_var, times)
end

function _to_lab_frame_velocity(soln, vars, times)
    timetrace = zeros(length(times))
    for var in vars
        val = real(substitute_all(unwrap(_remove_brackets(var)), soln))
        ω = real(real(unwrap(substitute_all(var.ω, soln))))
        if var.type == "u"
            timetrace .+= -ω * val * sin.(ω * times)
        elseif var.type == "v"
            timetrace .+= ω * val * cos.(ω * times)
        end
    end
    return timetrace
end

function _to_lab_frame(soln, vars, times)::Vector{AbstractFloat}
    timetrace = zeros(length(times))

    for var in vars
        val = real(substitute_all(unwrap(_remove_brackets(var)), soln))
        ω = real(unwrap(substitute_all(var.ω, soln)))
        if var.type == "u"
            timetrace .+= val * cos.(ω * times)
        elseif var.type == "v"
            timetrace .+= val * sin.(ω * times)
        elseif var.type == "a"
            timetrace .+= val
        end
    end
    return timetrace
end
