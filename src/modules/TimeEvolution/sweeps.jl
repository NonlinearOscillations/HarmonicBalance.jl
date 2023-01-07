export ParameterSweep

function ParameterSweep(functions::Dict, timespan::Tuple)
    t0, t1 = timespan[1], timespan[2]
    sweep_func=Dict{Num,Any}([])
    for swept_p in keys(functions)
        bounds = functions[swept_p]
        tfunc = swept_function(bounds, timespan)
        setindex!(sweep_func,tfunc,swept_p)
    end
    return ParameterSweep(sweep_func)
end


ParameterSweep(functions, timespan::Tuple) = ParameterSweep(Dict(functions), timespan)


function swept_function(bounds, timespan)
    t0, t1 = timespan
    function f(t)
        if t > t1
            bounds[2]
        elseif t < t0
            bounds[1]
        else
            ((t1 - t)/(t1-t0)) * bounds[1] + ((t - t0)/(t1-t0) * bounds[2])
        end
    end
    return f
end
