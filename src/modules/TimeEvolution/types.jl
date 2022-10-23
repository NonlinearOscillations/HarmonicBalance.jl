import Base: keys, getindex, +
export +, getindex

"""

Represents a sweep of one or more parameters of a `HarmonicEquation`.
During a sweep, the selected parameters vary linearly over some timespan and are constant elsewhere.

Sweeps of different variables can be combined using `+`.

# Fields
$(TYPEDFIELDS)

## Examples
```julia-repl
# create a sweep of parameter a from 0 to 1 over time 0 -> 100
julia> @variables a,b;
julia> sweep = ParameterSweep(a => [0., 1.], (0, 100));
julia> sweep[a](50)
0.5
julia> sweep[a](200)
1.0

# do the same, varying two parameters simultaneously
julia> sweep = ParameterSweep([a => [0.,1.], b => [0., 1.]], (0,100))
```
"""
struct ParameterSweep
    """Maps each swept parameter to a function."""
    functions::Dict{Num, Function}

    ParameterSweep(functions...) = new(Dict(functions...))
    ParameterSweep() = ParameterSweep([])

end


# overload so that ParameterSweep can be accessed like a Dict
keys(s::ParameterSweep) = keys(s.functions)
getindex(s::ParameterSweep, i) = getindex(s.functions, i)


# overload +
function +(s1::ParameterSweep, s2::ParameterSweep)
    common_params = intersect(keys(s1), keys(s2))
    !isempty(common_params) && error("cannot combine sweeps of the same parameter")
    return ParameterSweep(merge(s1.functions, s2.functions))

    # combine sweeps of the same parameter
    #interval_overlap(s1.timespan, s2.timespan) && error("cannot combine sweeps with overlapping timespans")
    #new_funcs = filter(x -> !in(x.first, common_params), all_params)
end


#interval_overlap((t1, t2), (t3, t4)) = t4 > t2 > t3 || t2 > t4 > t1