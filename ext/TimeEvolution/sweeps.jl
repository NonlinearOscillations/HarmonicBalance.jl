using HarmonicBalance: AdiabaticSweep

function HarmonicBalance.AdiabaticSweep(functions::Dict, timespan::Tuple)
    t0, t1 = timespan[1], timespan[2]
    sweep_func = Dict{Num,Any}([])
    for swept_p in keys(functions)
        bounds = functions[swept_p]
        tfunc = swept_function(bounds, timespan)
        setindex!(sweep_func, tfunc, swept_p)
    end
    return AdiabaticSweep(sweep_func)
end

function HarmonicBalance.AdiabaticSweep(functions, timespan::Tuple)
    return AdiabaticSweep(Dict(functions), timespan)
end

function swept_function(bounds, timespan)
    t0, t1 = timespan
    function f(t)
        if t > t1
            bounds[2]
        elseif t < t0
            bounds[1]
        else
            ((t1 - t) / (t1 - t0)) * bounds[1] + ((t - t0) / (t1 - t0) * bounds[2])
        end
    end
    return f
end

# overload so that AdiabaticSweep can be accessed like a Dict
Base.keys(s::AdiabaticSweep) = keys(s.functions)
Base.getindex(s::AdiabaticSweep, i) = getindex(s.functions, i)

# overload +
function Base.:+(s1::AdiabaticSweep, s2::AdiabaticSweep)
    common_params = intersect(keys(s1), keys(s2))
    !isempty(common_params) && error("cannot combine sweeps of the same parameter")
    return AdiabaticSweep(merge(s1.functions, s2.functions))

    # combine sweeps of the same parameter
    #interval_overlap(s1.timespan, s2.timespan) && error("cannot combine sweeps with overlapping timespans")
    #new_funcs = filter(x -> !in(x.first, common_params), all_params)
end
