const _set_Plots_default = Dict{Symbol,Any}([
    :fontfamily => "computer modern",
    :titlefont => "computer modern",
    :tickfont => "computer modern",
    :linewidth => 2,
    :legend_position => :outerright,
])

dim(res::Result) = length(size(res.solutions)) # give solution dimensionality

"""
$(TYPEDSIGNATURES)

## Plot a `Result` object.

Class selection done by passing `String` or `Vector{String}` as kwarg:

    class       :   only plot solutions in this class(es) ("all" --> plot everything)
    not_class   :   do not plot solutions in this class(es)

Other kwargs are passed onto Plots.gr().

See also `plot!`

The x,y,z arguments are Strings compatible with Symbolics.jl, e.g., `y=2*sqrt(u1^2+v1^2)` plots
the amplitude of the first quadratures multiplied by 2.

## 1D plots
    plot(res::Result; x::String, y::String, class="default", not_class=[], kwargs...)
    plot(res::Result, y::String; kwargs...) # take x automatically from Result

Default behaviour is to plot stable solutions as full lines, unstable as dashed.

If a sweep in two parameters were done, i.e., `dim(res)==2`, a one dimensional cut can be
plotted by using the keyword `cut` were it takes a `Pair{Num, Float64}` type entry. For example,
`plot(res, y="sqrt(u1^2+v1^2), cut=(λ => 0.2))` plots a cut at `λ = 0.2`.
###

## 2D plots
    plot(res::Result; z::String, branch::Int64, class="physical", not_class=[], kwargs...)

To make the 2d plot less chaotic it is required to specify the specific `branch` to plot, labeled by a `Int64`.

The x and y axes are taken automatically from `res`
"""
function Plots.plot(
    res::Result, varargs...; cut=Pair(missing, missing), kwargs...
)::Plots.Plot
    if dim(res) == 1
        plot1D(res, varargs...; _set_Plots_default..., kwargs...)
    elseif dim(res) == 2
        if ismissing(cut.first)
            plot2D(res, varargs...; _set_Plots_default..., kwargs...)
        else
            plot2D_cut(res, varargs...; cut=cut, _set_Plots_default..., kwargs...)
        end
    else
        error("Data dimension ", dim(res), " not supported")
    end
end

"""
$(TYPEDSIGNATURES)

Similar to `plot` but adds a plot onto an existing plot.
"""
function Plots.plot!(res::Result, varargs...; kwargs...)::Plots.Plot
    return plot(res, varargs...; add=true, _set_Plots_default..., kwargs...)
end
"""
$(TYPEDSIGNATURES)

Return an array of bools to mark solutions in `res` which fall into `classes` but not `not_classes`.
Only `branches` are considered.
"""
function _get_mask(res, classes, not_classes=[]; branches=1:branch_count(res))
    classes == "all" && return fill(trues(length(branches)), size(res.solutions))
    bools = vcat(
        [res.classes[c] for c in _str_to_vec(classes)],
        [map(.!, res.classes[c]) for c in _str_to_vec(not_classes)],
    )
    #m = map( x -> [getindex(x, b) for b in [branches...]], map(.*, bools...))

    return m = map(x -> x[[branches...]], map(.*, bools...))
end

"""
$(TYPEDSIGNATURES)

Go over a solution and an equally-sized array (a "mask") of booleans.
true  -> solution unchanged
false -> changed to NaN (omitted from plotting)
"""
function _apply_mask(solns::Array{Vector{ComplexF64}}, booleans)
    factors = replace.(booleans, 0 => NaN)
    return map(.*, solns, factors)
end
function _apply_mask(solns::Vector{Vector{Vector{ComplexF64}}}, booleans)
    Nan_vector = NaN .* similar(solns[1][1])
    new_solns = [
        [booleans[i][j] ? solns[i][j] : Nan_vector for j in eachindex(solns[i])] for
        i in eachindex(solns)
    ]
    return new_solns
end

""" Project the array `a` into the real axis, warning if its contents are complex. """
function _realify(a::Array{T} where {T<:Number}; warning="")
    warned = false
    a_real = similar(a, Float64)
    for i in eachindex(a)
        if !isnan(a[i]) && !warned && !is_real(a[i])
            @warn "Values with non-negligible complex parts have
            been projected on the real axis! " * warning
            warned = true
        end
        a_real[i] = real(a[i])
    end
    return a_real
end

_str_to_vec(s::Vector) = s
_str_to_vec(s) = [s]

# return true if p already has a label for branch index idx
function _is_labeled(p::Plot, idx::Int64)
    return in(string(idx), [sub[:label] for sub in p.series_list])
end

function plot1D(
    res::Result;
    x::String="default",
    y::String,
    class="default",
    not_class=[],
    branches=1:branch_count(res),
    add=false,
    kwargs...,
)::Plots.Plot
    if class == "default"
        args = [:x => x, :y => y, :branches => branches]
        if not_class == [] # plot stable full, unstable dashed
            p = plot1D(res; args..., class=["physical", "stable"], add=add, kwargs...)
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
            p = plot1D(
                res; args..., not_class=not_class, class="physical", add=add, kwargs...
            )
            return p
        end
    end

    dim(res) != 1 &&
        error("The results are two dimensional. Consider using the `cut` keyword.")
    x = x == "default" ? string(first(keys(res.swept_parameters))) : x
    X = transform_solutions(res, x; branches=branches)
    Y = transform_solutions(res, y; branches=branches, realify=true)
    Y = _apply_mask(Y, _get_mask(res, class, not_class; branches=branches))

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

plot1D(res::Result, y::String; kwargs...) = plot1D(res; y=y, kwargs...)

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
    Z =
        getindex.(
            _apply_mask(
                transform_solutions(res, z; branches=branch, realify=true),
                _get_mask(res, class, not_class; branches=branch),
            ),
            1,
        ) # there is only one branch
    p = add ? Plots.plot!() : Plots.plot() # start a new plot if needed

    ylab, xlab = latexify.(string.(keys(res.swept_parameters)))
    return p = plot!(
        map(_realify, [Float64.(Y), Float64.(X), Z])...;
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
            p = plot2D_cut(
                res; y=y, cut=cut, class=["physical", "stable"], add=add, kwargs...
            )
            plot2D_cut(
                res;
                y=y,
                cut=cut,
                class="physical",
                not_class="stable",
                add=true,
                style=:dash,
                kwargs...,
            )
            return p
        else
            p = plot2D_cut(
                res; y=y, cut=cut, not_class=not_class, class="physical", add=add, kwargs...
            )
            return p
        end
    end

    cut_par, cut_value = cut
    # compare strings beacuse type Num cannot be compared
    swept_pars = res.swept_parameters.keys
    x_index = findfirst(sym -> string(sym) != string(cut_par), swept_pars)
    isnothing(x_index) && error("The variable $cut_par was not swept over.")
    x = swept_pars[x_index]

    # the swept params are ranges and thus a sorted search can be performed
    cut_par_index = searchsortedfirst(res.swept_parameters[cut_par], cut_value)
    if !(cut_par_index ∈ 1:length(res.swept_parameters[cut_par]))
        throw(
            ArgumentError(
                "The value $cut_value is not found in the swept range of $cut_par."
            ),
        )
    end

    X = res.swept_parameters[x]
    Y = _apply_mask(
        transform_solutions(res, y; realify=true), _get_mask(res, class, not_class)
    ) # first transform, then filter
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
            Float64.(X),
            _realify(getindex.(branches, k));
            color=k,
            label=l,
            xlabel=latexify(string(x)),
            ylabel=latexify(y),
            kwargs...,
        )
    end

    return p
end

plot2D(res::Result, z::String; kwargs...) = plot2D(res; z=z, kwargs...)

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
function plot_phase_diagram(res::Result; kwargs...)::Plots.Plot
    if dim(res) == 1
        plot_phase_diagram_1D(res; _set_Plots_default..., kwargs...)
    elseif dim(res) == 2
        plot_phase_diagram_2D(res; _set_Plots_default..., kwargs...)
    else
        error("Data dimension ", dim(res), " not supported")
    end
end

function plot_phase_diagram(res::Result, class::String; kwargs...)
    return plot_phase_diagram(res; class=class, kwargs...)
end

function plot_phase_diagram_2D(
    res::Result; class="physical", not_class=[], kwargs...
)::Plots.Plot
    X, Y = values(res.swept_parameters)
    Z = sum.(_get_mask(res, class, not_class))

    xlab, ylab = latexify.(string.(keys(res.swept_parameters)))

    # cannot set heatmap ticks (Plots issue #3560)
    return heatmap(
        Float64.(X),
        Float64.(Y),
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
    X = first(values(res.swept_parameters))
    Y = sum.(_get_mask(res, class, not_class))
    return plot(
        Float64.(X),
        Y;
        xlabel=latexify(string(keys(res.swept_parameters)...)),
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
    plot_spaghetti(res::Result; x, y, z, kwargs...)

Plot a three dimension line plot of a `Result` object as a function of the parameters.
Works with 1D and 2D datasets.

Class selection done by passing `String` or `Vector{String}` as kwarg:

    class::String       :   only count solutions in this class ("all" --> plot everything)
    not_class::String   :   do not count solutions in this class

Other kwargs are passed onto Plots.gr()
"""
function plot_spaghetti(
    res::Result;
    x::String,
    y::String,
    z::String,
    class="default",
    not_class=[],
    add=false,
    kwargs...,
)::Plots.Plot
    if dim(res) == 2
        error("Data dimension ", dim(res), " not supported")
    end

    if class == "default"
        if not_class == [] # plot stable full, unstable dashed
            p = plot_spaghetti(
                res; x=x, y=y, z=z, class=["physical", "stable"], add=add, kwargs...
            )
            plot_spaghetti(
                res;
                x=x,
                y=y,
                z=z,
                class="physical",
                not_class="stable",
                add=true,
                style=:dash,
                kwargs...,
            )
            return p
        else
            p = plot_spaghetti(
                res;
                x=x,
                y=y,
                z=z,
                class="physical",
                not_class=not_class,
                add=add,
                kwargs...,
            )
            return p
        end
    end

    vars = res.problem.variables
    x_index = findfirst(sym -> string(sym) == x, vars)
    y_index = findfirst(sym -> string(sym) == y, vars)
    isnothing(x_index) && error("The variable $x is not a defined variable.")
    isnothing(y_index) && error("The variable $y is not a defined variable.")

    swept_pars = res.swept_parameters.keys
    z_index = findfirst(sym -> string(sym) == z, swept_pars)
    isnothing(z_index) && error("The variable $z was not swept over.")

    Z = res.swept_parameters.vals[z_index]
    X = _apply_mask(transform_solutions(res, x), _get_mask(res, class, not_class))
    Y = _apply_mask(transform_solutions(res, y), _get_mask(res, class, not_class))

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
            Float64.(Z);
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
