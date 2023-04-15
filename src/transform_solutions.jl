export transform_solutions


_parse_expression(exp) = exp isa String ? Num(eval(Meta.parse(exp))) : exp


"""
$(TYPEDSIGNATURES)

Takes a `Result` object and a string `f` representing a Symbolics.jl expression.
Returns an array with the values of `f` evaluated for the respective solutions.
Additional substitution rules can be specified in `rules` in the format `("a" => val)` or `(a => val)`
"""
function transform_solutions(res::Result, f::String; branches = 1:branch_count(res), rules=Dict(), target_type=ComplexF64)
    # a string is used as input - a macro would not "see" the user's namespace while the user's namespace does not "see" the variables
    transformed = [Vector{target_type}(undef, length(branches)) for k in res.solutions] # preallocate

    func = _build_substituted(res, f; rules=rules)

    # preallocate an array for the numerical values, rewrite parts of it
    # when looping through the solutions
    n_vars = length(get_variables(res))
    n_pars = length(res.swept_parameters)
    vals = Vector{ComplexF64}(undef, n_vars + n_pars)

    for idx in CartesianIndices(res.solutions)
        params_values = res.swept_parameters[Tuple(idx)...]
        vals[end-n_pars+1:end] .= params_values # param values are common to all branches
        for (k, branch) in enumerate(branches)
            vals[1:n_vars] .= res.solutions[idx][branch]
            transformed[idx][k] = Base.invokelatest(func, vals)
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


""" Parse `expr` into a Symbolics.jl expression, substitute with `rules` and build a function taking free_symbols """
function _build_substituted(expr::String, rules, free_symbols)

    subbed = substitute_all(_parse_expression(expr), rules)
    comp_func = build_function(subbed, free_symbols)

    return eval(comp_func)
end

""" Parse `expr` into a Symbolics.jl expression, substituting the fixed parameters of `res`
The resulting function takes in the values of the variables and swept parameters. """
function _build_substituted(res::Result, expr; rules=Dict())

   # define variables in rules in this namespace
   new_keys = declare_variable.(string.(keys(Dict(rules))))
   fixed_subs = merge(res.fixed_parameters, Dict(zip(new_keys, values(Dict(rules)))))

   free_symbols = vcat(res.problem.variables, collect(keys(res.swept_parameters)))
   return _build_substituted(expr, fixed_subs, free_symbols)

end


## move masks here
