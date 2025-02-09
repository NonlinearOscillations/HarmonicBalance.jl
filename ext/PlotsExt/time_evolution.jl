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
    soln::SciMLBase.ODESolution,
    funcs,
    harm_eq::HarmonicBalance.HarmonicEquation;
    add=false,
    kwargs...,
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

function Plots.plot!(soln::SciMLBase.ODESolution, varargs...; kwargs...)
    return plot(soln, varargs...; add=true, kwargs...)
end

"""
$(TYPEDSIGNATURES)

Plot a bifurcation diagram from a continuation sweep starting from `starting_branch` using
the [Result](@ref) struct `res`. Time integration is used to determine what follow up branch
in the continuation is.

# Keyword arguments
- `x::String`: Expression for the x-axis variable
- `y::String`: Expression for the y-axis variable
- `sweep::String="right"`: Direction to follow the branch ("right" or "left")
- `tf::Real=10000`: Final time for time integration
- `ϵ::Real=1e-4`: Tolerance for branch following
- `kwargs...`: Additional plotting arguments passed to Plots.jl
- Class selection done by passing `String` or `Vector{String}` as kwarg:

    class::String       :   only count solutions in this class ("all" --> plot everything)
    not_class::String   :   do not count solutions in this class


# Returns
- A Plots.jl plot object containing the bifurcation diagram with the followed branch

# Description
This function creates a bifurcation diagram using [`follow_branch`](@ref).
The followed branch is plotted as a dashed gray line.
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
    p = plot(res; x, y, class, not_class, kwargs...)

    followed_branch, Ys = HarmonicBalance.follow_branch(
        starting_branch, res; y, sweep=sweep, tf=tf, ϵ=ϵ
    )
    Y_followed = [
        Ys[param_idx][branch] for (param_idx, branch) in enumerate(followed_branch)
    ]
    X = swept_parameter(res, x)

    Plots.plot!(
        p,
        real.(X),
        real.(Y_followed);
        linestyle=:dash,
        c=:gray,
        label=sweep * " sweep",
        _set_Plots_default...,
        kwargs...,
    )
    return p
end
