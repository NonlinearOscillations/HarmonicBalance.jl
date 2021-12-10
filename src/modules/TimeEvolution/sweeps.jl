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


#=

"""
Returns functions corresponding to multiple linear sweeps. sweeps should contain a dictionary with the swept variable and its initial
intermediate and final values. timespan contains the initial, intermediate and final time values for the sweep.
"""
function multiple_sweeps_functions(sweeps, timespans)
    swept = keys(sweeps)
    reshaped_sweeps=reshape_sweeps(sweeps,timespans)
    nullfunction(t)=0
    sweep_funcs = Dict{Num,Any}(zip(swept,[nullfunction for i in 1:length(swept)]))
    for i in 1:length(timespans)-1
        temp_sweep = sweep_functions(reshaped_sweeps[i], timespans[i:i+1])
        for swept_p in swept
            sweep_funcs[swept_p]= let temp = sweep_funcs[swept_p];t -> temp(t)+interval(t,timespans[i],timespans[i+1])*temp_sweep[swept_p](t); end
        end
    end
    return sweep_funcs
end


function heaviside0(x)
    heaviside = 0
    if (x > 0) heaviside = 1
    end
    return heaviside
end
function heaviside1(x)
    heaviside = 0
    if (x >= 0) heaviside = 1
    end
    return heaviside
end

function interval(t, a, b)
   heaviside1(t-a) - heaviside1(t-b)
end


"""
Seperates the sweeps provided in sweeps into an array of separate sweeps.
E.g. Dict(ω=>[0.95;1.;1.;3],λ => [5E-2;5E-2;1E-2]) is transformed into [Dict(λ => [0.05, 0.05], ω => [0.95, 1.0]),
Dict(λ => [0.05, 0.01], ω => [1.0, 1.0])]
"""
function reshape_sweeps(sweeps, timespans)

    # error alert??
    swept = keys(sweeps)
    sweeps_reshaped=[Dict(zip(swept, [zeros(2) for i in  1:length(swept)])) for i in 1:length(timespans)-1]

    for swept_p in swept
        lines = sweeps[swept_p]
        for i in 1:length(timespans)-1
            sweeps_reshaped[i][swept_p]=lines[i:i+1]
        end
    end
    return sweeps_reshaped
end




function sweep_functions(sweep, timespan)
    " Returns functions corresponding to a linear sweep. sweep should contain a dictionary with the swept variable and its initial
      and final value. timespan contains the initial and final time value for the sweep."
    t0, t1 = timespan[1], timespan[2]
    swept = keys(sweep)
    sweep_func=Dict{Num,Any}([])
    for swept_p in swept
        line = sweep[swept_p]
        tfunc(T) = ((t1 - T)/(t1-t0)) * line[1] + ((T - t0)/(t1-t0) * line[2])
        setindex!(sweep_func,tfunc,swept_p)
    end
    return sweep_func
end
=#