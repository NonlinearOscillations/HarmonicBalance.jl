
"""
$(TYPEDSIGNATURES)

## Plot a `Result` object.

Class selection done by passing `String` or `Vector{String}` as kwarg:

    class       :   only plot solutions in this class(es) ("all" --> plot everything)
    not_class   :   do not plot solutions in this class(es)

Other kwargs are passed onto Plots.gr().

See also `plot!`

The x,y,z arguments are Strings compatible with Symbolics.jl, e.g., `y=2*sqrt(u1^2+v1^2)`
plots the amplitude of the first quadratures multiplied by 2.

## 1D plots
    plot(res::Result; x::String, y::String, class="default", not_class=[], kwargs...)
    plot(res::Result, y::String; kwargs...) # take x automatically from Result

Default behaviour is to plot stable solutions as full lines, unstable as dashed.

If a sweep in two parameters were done, i.e., `dimension(res)==2`, a one dimensional cut can
be plotted by using the keyword `cut` were it takes a `Pair{Num, Float}` type entry.
For example, `plot(res, y="sqrt(u1^2+v1^2), cut=(λ => 0.2))` plots a cut at `λ = 0.2`.

## 2D plots
    plot(res::Result; z::String, branch::Int64, class="physical", not_class=[], kwargs...)

To make the 2d plot less chaotic it is required to specify the specific `branch` to plot,
labeled by a `Int64`.

The x and y axes are taken automatically from `res`
"""
function Plots.plot(
    res::Result{D,S,P}, varargs...; cut=Pair(missing, missing), kwargs...
)::Plots.Plot where {D,S,P}
    if D == 1
        plot1D(res, varargs...; _set_Plots_default..., kwargs...)
    elseif D == 2
        if ismissing(cut.first)
            plot2D(res, varargs...; _set_Plots_default..., kwargs...)
        else
            plot2D_cut(res, varargs...; cut=cut, _set_Plots_default..., kwargs...)
        end
    else
        error("Data dimension ", D, " not supported")
    end
end

"""
$(TYPEDSIGNATURES)

Similar to `plot` but adds a plot onto an existing plot.
"""
function Plots.plot!(res::Result, varargs...; kwargs...)::Plots.Plot
    return plot(res, varargs...; add=true, _set_Plots_default..., kwargs...)
end

# return true if p already has a label for branch index idx
function _is_labeled(p::Plots.Plot, idx::Int64)
    return in(string(idx), [sub[:label] for sub in p.series_list])
end

function plot1D(
    res::Result{D};
    x::String="default",
    y::String,
    class="default",
    not_class=[],
    branches=1:branch_count(res),
    add=false,
    kwargs...,
)::Plots.Plot where {D}
    if class == "default"
        args = [:x => x, :y => y, :branches => branches]
        if not_class == [] # plot stable full, unstable dashed
            p = plot1D(res; args..., class=["physical", "stable"], add, kwargs...)
            plot1D(
                res;
                args...,
                class="physical",
                not_class="stable",
                add=true,
                style=:dash,
                kwargs...,
            )
            return p
        else
            p = plot1D(res; args..., not_class, class="physical", add, kwargs...)
            return p
        end
    end

    D != 1 && error("The results are two dimensional. Consider using the `cut` keyword.")
    x = x == "default" ? string(first(keys(res.swept_parameters))) : x
    X = transform_solutions(res, x; branches)
    Y = get_solutions(res, y; branches, realify=true, class, not_class)

    # reformat and project onto real, warning if needed
    branch_data = [
        _realify(getindex.(Y, i); warning="branch " * string(k)) for
        (i, k) in enumerate(branches)
    ]

    # start a new plot if needed
    p = add ? Plots.plot!() : Plots.plot()

    # colouring is matched to branch index - matched across plots
    for k in findall(x -> !all(isnan.(x)), branch_data) # skip NaN branch_data
        global_index = branches[k]
        lab = _is_labeled(p, global_index) ? nothing : global_index
        Plots.plot!(
            _realify(getindex.(X, k)),
            branch_data[k];
            color=k,
            label=lab,
            xlabel=latexify(x),
            ylabel=latexify(y),
            kwargs...,
        )
    end

    return p
end

plot1D(res::Result, y::String; kwargs...) = plot1D(res; y, kwargs...)

function plot2D(
    res::Result;
    z::String,
    branch::Int64,
    class="physical",
    not_class=[],
    add=false,
    kwargs...,
)::Plots.Plot
    X, Y = values(res.swept_parameters)
    Z = getindex.(get_solutions(res, z; branches=branch, realify=true, class, not_class), 1) # there is only one branch
    p = add ? Plots.plot!() : Plots.plot() # start a new plot if needed

    ylab, xlab = get_labels(res)
    return p = plot!(
        map(_realify, [real.(Y), real.(X), Z])...;
        st=:surface,
        color=:blue,
        opacity=0.5,
        xlabel=xlab,
        ylabel=ylab,
        zlabel=latexify(z),
        colorbar=false,
        kwargs...,
    )
end

function plot2D_cut(
    res::Result; y::String, cut::Pair, class="default", not_class=[], add=false, kwargs...
)
    if class == "default"
        if not_class == [] # plot stable full, unstable dashed
            p = plot2D_cut(res; y, cut=cut, class=["physical", "stable"], add, kwargs...)
            plot2D_cut(
                res;
                y,
                cut=cut,
                class="physical",
                not_class="stable",
                add=true,
                style=:dash,
                kwargs...,
            )
            return p
        else
            p = plot2D_cut(res; y, cut=cut, not_class, class="physical", add, kwargs...)
            return p
        end
    end

    cut_par, cut_value = cut
    is_swept_parameter(res, cut_par)

    # the swept params are ranges and thus a sorted search can be performed
    cut_par_index = searchsortedfirst(res.swept_parameters[cut_par], cut_value)
    if !(cut_par_index ∈ 1:length(res.swept_parameters[cut_par]))
        throw(
            ArgumentError(
                "The value $cut_value is not found in the swept range of $cut_par."
            ),
        )
    end

    X = swept_parameter(res, cut_par)
    Y = get_solutions(res, y; realify=true, class, not_class) # first transform, then filter

    x_index = findfirst(sym -> string(sym) == string(cut_par), res.swept_parameters.keys)
    branches = x_index == 1 ? Y[:, cut_par_index] : Y[cut_par_index, :]

    branch_data = [
        _realify(getindex.(branches, i); warning="branch " * string(k)) for
        (i, k) in enumerate(1:branch_count(res))
    ]

    # start a new plot if needed
    p = add ? plot!() : plot()

    # colouring is matched to branch index - matched across plots
    for k in findall(branch -> !all(isnan.(branch)), branch_data) # skip NaN branches but keep indices
        l = _is_labeled(p, k) ? nothing : k
        Plots.plot!(
            real.(X),
            _realify(getindex.(branches, k));
            color=k,
            label=l,
            xlabel=latexify(string(cut_par)),
            ylabel=latexify(y),
            kwargs...,
        )
    end

    return p
end

plot2D(res::Result, z::String; kwargs...) = plot2D(res; z, kwargs...)

###
# PHASE DIAGRAMS
###

"""
$(TYPEDSIGNATURES)

Plot the number of solutions in a `Result` object as a function of the parameters.
Works with 1D and 2D datasets.

Class selection done by passing `String` or `Vector{String}` as kwarg:

    class::String       :   only count solutions in this class ("all" --> plot everything)
    not_class::String   :   do not count solutions in this class

Other kwargs are passed onto Plots.gr()
"""
function HarmonicBalance.plot_phase_diagram(res::Result{D}; kwargs...)::Plots.Plot where {D}
    if D == 1
        plot_phase_diagram_1D(res; _set_Plots_default..., kwargs...)
    elseif D == 2
        plot_phase_diagram_2D(res; _set_Plots_default..., kwargs...)
    else
        error("Data dimension ", D, " not supported")
    end
end

function HarmonicBalance.plot_phase_diagram(res::Result, class::String; kwargs...)
    return plot_phase_diagram(res; class, kwargs...)
end

function plot_phase_diagram_2D(
    res::Result; class="physical", not_class=[], kwargs...
)::Plots.Plot
    X, Y = swept_parameters(res)
    Z = phase_diagram(res; class, not_class)

    xlab, ylab = get_labels(res)

    # cannot set heatmap ticks (Plots issue #3560)
    return heatmap(
        real.(X),
        real.(Y),
        transpose(Z);
        xlabel=xlab,
        ylabel=ylab,
        color=:viridis,
        kwargs...,
    )
end

function plot_phase_diagram_1D(
    res::Result; class="physical", not_class=[], kwargs...
)::Plots.Plot
    X = swept_parameters(res)
    Y = phase_diagram(res; class, not_class)
    return plot(
        X,
        Y;
        xlabel=get_labels(res),
        ylabel="#",
        legend=false,
        yticks=1:maximum(Y),
        kwargs...,
    )
end

###
# Spaghetti Plot
###

"""
$(TYPEDSIGNATURES)

Plot a three dimension line plot of a `Result` object as a function of the parameters.
Works with 1D and 2D datasets.

Class selection done by passing `String` or `Vector{String}` as kwarg:

    class::String       :   only count solutions in this class ("all" --> plot everything)
    not_class::String   :   do not count solutions in this class

Other kwargs are passed onto Plots.gr()
"""
function HarmonicBalance.plot_spaghetti(
    res::Result{D};
    x::String,
    y::String,
    z::String,
    class="default",
    not_class=[],
    add=false,
    kwargs...,
)::Plots.Plot where {D}
    if D == 2
        error("Data dimension ", D, " not supported")
    end

    if class == "default"
        if not_class == [] # plot stable full, unstable dashed
            p = HarmonicBalance.plot_spaghetti(
                res; x, y, z, class=["physical", "stable"], add, kwargs...
            )
            HarmonicBalance.plot_spaghetti(
                res;
                x,
                y,
                z,
                class="physical",
                not_class="stable",
                add=true,
                style=:dash,
                kwargs...,
            )
            return p
        else
            p = HarmonicBalance.plot_spaghetti(
                res; x, y, z, class="physical", not_class, add, kwargs...
            )
            return p
        end
    end

    is_variable(res, x)
    is_variable(res, y)

    Z = swept_parameter(res, z)
    X = get_solutions(res, x; class, not_class)
    Y = get_solutions(res, y; class, not_class)

    # start a new plot if needed
    p = add ? Plots.plot!() : Plots.plot()

    branch_data = [
        _realify(getindex.(X, i); warning="branch " * string(k)) for
        (i, k) in enumerate(1:branch_count(res))
    ]

    # colouring is matched to branch index - matched across plots
    for k in findall(x -> !all(isnan.(x)), branch_data) # skip NaN branches but keep indices
        l = _is_labeled(p, k) ? nothing : k
        Plots.plot!(
            _realify(getindex.(X, k)),
            _realify(getindex.(Y, k)),
            real.(Z);
            _set_Plots_default...,
            color=k,
            label=l,
            xlabel=latexify(x),
            ylabel=latexify(y),
            zlabel=latexify(z),
            xlim=:symmetric,
            ylim=:symmetric,
            kwargs...,
        )
    end
    return p
end
