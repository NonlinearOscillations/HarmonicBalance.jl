using Plots: Plots, plot, plot!
using Latexify: latexify

"""
    plot(soln::ODESolution, f::String, harm_eq::HarmonicEquation; kwargs...)

Plot a function `f` of a time-dependent solution `soln` of `harm_eq`.

## As a function of time

    plot(soln::ODESolution, f::String, harm_eq::HarmonicEquation; kwargs...)

`f` is parsed by Symbolics.jl

## parametric plots

    plot(soln::ODESolution, f::Vector{String}, harm_eq::HarmonicEquation; kwargs...)

Parametric plot of f[1] against f[2]

Also callable as plot!
"""
function Plots.plot(
    soln::OrdinaryDiffEq.ODESolution, funcs, harm_eq::HarmonicEquation; add=false, kwargs...
)

    # start a new plot if needed
    p = add ? plot!() : plot()

    if funcs isa String || length(funcs) == 1
        plot!(
            soln.t,
            transform_solutions(soln, funcs, harm_eq);
            _set_Plots_default...,
            xlabel="time",
            ylabel=latexify(funcs),
            legend=false,
            kwargs...,
        )
    elseif length(funcs) == 2 # plot of func vs func
        plot!(
            transform_solutions(soln, funcs, harm_eq)...;
            _set_Plots_default...,
            xlabel=latexify(funcs[1]),
            ylabel=latexify(funcs[2]),
            legend=false,
            kwargs...,
        )
    else
        error("Invalid plotting argument: ", funcs)
    end
end

function Plots.plot!(soln::OrdinaryDiffEq.ODESolution, varargs...; kwargs...)
    return plot(soln, varargs...; add=true, kwargs...)
end

"""
1D plot with the followed branch highlighted
"""
function HarmonicBalance.plot_1D_solutions_branch(
    starting_branch::Int64,
    res::Result;
    x::String,
    y::String,
    sweep="right",
    tf=10000,
    ϵ=1e-4,
    class="default",
    not_class=[],
    kwargs...,
)
    p = plot(res; x=x, y=y, class=class, not_class=not_class, kwargs...)

    followed_branch, Ys = follow_branch(starting_branch, res; y=y, sweep=sweep, tf=tf, ϵ=ϵ)
    Y_followed = [
        Ys[param_idx][branch] for (param_idx, branch) in enumerate(followed_branch)
    ]
    X = real.(res.swept_parameters[_parse_expression(x)])

    Plots.plot!(
        p,
        X,
        real.(Y_followed);
        linestyle=:dash,
        c=:gray,
        label=sweep * " sweep",
        _set_Plots_default...,
        kwargs...,
    )
    return p
end
