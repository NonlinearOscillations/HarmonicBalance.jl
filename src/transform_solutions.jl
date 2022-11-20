export transform_solutions


_parse_expression(exp) = exp isa String ? Num(eval(Meta.parse(exp))) : exp


"""
$(TYPEDSIGNATURES)

Takes a `Result` object and a string `f` representing a Symbolics.jl expression.
Returns an array with the values of `f` evaluated for the respective solutions.
Additional substitution rules can be specified in `rules` in the format `("a" => val)` or `(a => val)`
"""
function transform_solutions(res::Result, f::String, branches = 1:branch_count(res); rules=Dict())
    # a string is used as input - a macro would not "see" the user's namespace while the user's namespace does not "see" the variables
    transformed = [Vector{ComplexF64}(undef, length(branches)) for k in res.solutions] # preallocate

    # define variables in rules in this namespace
    new_keys = declare_variable.(string.(keys(Dict(rules))))
    expr = f isa String ? _parse_expression(f) : f

    fixed_subs = merge(res.fixed_parameters, Dict(zip(new_keys, values(Dict(rules)))))
    expr = substitute_all(expr, Dict(fixed_subs))

    vars = res.problem.variables
    all_symbols = cat(vars, collect(keys(res.swept_parameters)), dims=1)
    comp_func = build_function(expr, all_symbols)
    f = eval(comp_func)

    # preallocate an array for the numerical values, rewrite parts of it
    # when looping through the solutions
    vals = Vector{ComplexF64}(undef, length(all_symbols))
    n_vars = length(vars)
    n_pars = length(all_symbols) - n_vars

    for idx in CartesianIndices(res.solutions)
        params_values = res.swept_parameters[Tuple(idx)...]
        vals[end-n_pars+1:end] .= params_values # param values are common to all branches
        for (k, branch) in enumerate(branches)
            vals[1:n_vars] .= res.solutions[idx][branch]
            transformed[idx][k] = Base.invokelatest(f, vals)
        end
    end
    return transformed
end

transform_solutions(res::Result, fs::Vector{String}; kwargs...) = [transform_solutions(res, f; kwargs...) for f in fs]

# a simplified version meant to work with arrays of solutions
# cannot parse parameter values -- meant for time-dependent results
function transform_solutions(soln::Vector, f::String, harm_eq::HarmonicEquation)

    vars = _remove_brackets(get_variables(harm_eq))
    transformed = Vector{ComplexF64}(undef, length(soln))

    # parse the input with Symbolics
    expr = HarmonicBalance._parse_expression(f)

    rule(u) = Dict(zip(vars, u))

    transformed = map( x -> substitute_all(expr, rule(x)), soln)
    return convert(typeof(soln[1]), transformed)
end