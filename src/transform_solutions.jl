export transform_solutions


_parse_expression(exp) = exp isa String ? Num(eval(Meta.parse(exp))) : exp


"""
$(TYPEDSIGNATURES)

Takes a `Result` object and a string `f` representing a Symbolics.jl expression.
Returns an array with the values of `f` evaluated for the respective solutions.
Additional substitution rules can be specified in `rules` in the format `("a" => val)` or `(a => val)`
"""
function transform_solutions(res::Result, func; branches = 1:branch_count(res), kwargs...)
    # preallocate an array for the numerical values, rewrite parts of it
    # when looping through the solutions
    pars = res.swept_parameters |> values |> collect
    n_vars = length(get_variables(res))
    n_pars = length(pars)

    vtype = isa(Base.invokelatest(func, rand(ComplexF64, n_vars+n_pars); kwargs...), Bool) ? BitVector : Vector{ComplexF64}
    transformed = _similar(vtype, res; branches=branches)

    batches = Iterators.partition(CartesianIndices(res.solutions), ceil(Int, length(res.solutions)/Threads.nthreads()))
    Threads.@threads for batch in batches |> collect
        _vals = Vector{ComplexF64}(undef, n_vars + n_pars)
        for idx in batch
            for i in 1:length(idx) # param values are common to all branches
                _vals[end-n_pars+i] = pars[i][idx[i]]
            end
            for (k, branch) in enumerate(branches)
                _vals[1:n_vars] .= res.solutions[idx][branch]
                transformed[idx][k] = Base.invokelatest(func, _vals; kwargs...) # beware, func may be mutating
            end
        end
    end
    return transformed
end

function transform_solutions(res::Result, f::String; rules=Dict(), kwargs...)
    # a string is used as input
    # a macro would not "see" the user's namespace while the user's namespace does not "see" the variables
    func = _build_substituted(f, res; rules=rules)
    transform_solutions(res, func; kwargs...)
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
function _build_substituted(expr, res::Result; rules=Dict())

   # define variables in rules in this namespace
   new_keys = declare_variable.(string.(keys(Dict(rules))))
   fixed_subs = merge(res.fixed_parameters, Dict(zip(new_keys, values(Dict(rules)))))

   return _build_substituted(expr, fixed_subs, _free_symbols(res))

end

function _similar(type, res::Result; branches=1:branch_count(res))
    [type(undef, length(branches)) for k in res.solutions]
end

## move masks here

###
# TRANSFORMATIONS TO THE LAB frame
###

function to_lab_frame(soln, res, times)

    timetrace = zeros(length(times))

    for var in res.problem.eom.variables
        val = real(substitute_all(_remove_brackets(var), soln))
        ω = substitute_all(var.ω, soln) |> Float64
        if var.type == "u"
            timetrace .+= val*cos.(ω * times)
        elseif var.type == "v"
            timetrace .+= val*sin.(ω * times)
        elseif var.type == "a"
            timetrace .+= val
        end
    end
    timetrace
end


"""
    to_lab_frame(res::Result, times; index::Int, branch::Int)
    to_lab_frame(soln::OrderedDict, res::Result, times)

Transform a solution into the lab frame (i.e., invert the harmonic ansatz) for `times`.
Either extract the solution from `res::Result` by `index` and `branch` or input `soln::OrderedDict` explicitly.
"""
to_lab_frame(res::Result, times; index::Int, branch::Int) = to_lab_frame(res[index][branch], res, times)


function to_lab_frame_velocity(soln, res, times)

    timetrace = zeros(length(times))

    for var in res.problem.eom.variables
        val = real(substitute_all(_remove_brackets(var), soln))
        ω = real(substitute_all(var.ω, soln)) |> Float64
        if var.type == "u"
            timetrace .+= -ω*val*sin.(ω * times)
        elseif var.type == "v"
            timetrace .+= ω*val*cos.(ω * times)
        end
    end
    timetrace
end


"""
    to_lab_frame_velocity(res::Result, times; index::Int, branch::Int)
    to_lab_frame_velocity(soln::OrderedDict, res::Result, times)

Transform a solution's velocity into the lab frame (i.e., invert the harmonic ansatz for dx/dt ) for `times`.
Either extract the solution from `res::Result` by `index` and `branch` or input `soln::OrderedDict` explicitly.
"""
to_lab_frame_velocity(res::Result, times; index, branch) = to_lab_frame_velocity(res[index][branch], res, times)
