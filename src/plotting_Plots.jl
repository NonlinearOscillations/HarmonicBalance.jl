using Plots, Latexify
import Plots.plot, Plots.plot!; export plot, plot!, plot_phase_diagram, savefig

const _set_Plots_default = Dict{Symbol, Any}([
    :fontfamily => "computer modern",
    :titlefont => "computer modern",
    :tickfont => "computer modern",
    :linewidth => 2,
    :legend => :outerright])



dim(res::Result) = length(size(res.solutions)) # give solution dimensionality


"""
$(TYPEDSIGNATURES)

## Plot a `Result` object.

Class selection done by passing `String` or `Vector{String}` as kwarg:

    class       :   only plot solutions in this class(es) ("all" --> plot everything)
    not_class   :   do not plot solutions in this class(es)

Other kwargs are passed onto Plots.gr()

See also `plot!`

The x,y,z arguments are Strings compatible with Symbolics.jl

## 1D plots
    plot(res::Result; x::String, y::String, class="default", not_class=[], kwargs...)
    plot(res::Result, y::String; kwargs...) # take x automatically from Result

Default behaviour is to plot stable solutions as full lines, unstable as dashed
###

## 2D plots
    plot(res::Result; z::String, branch::Int64, class="physical", not_class=[], kwargs...)

To make the 2d plot less chaotic it is required to specify the specific `branch` to plot, labeled by a `Int64`.

The x and y axes are taken automatically from `res`
"""
function plot(res::Result, varargs...; cut=Pair(missing, missing), kwargs...)::Plots.Plot
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
function plot!(res::Result, varargs...; kwargs...)::Plots.Plot
    plot(res, varargs...; add=true, _set_Plots_default..., kwargs...)
end
"""
$(TYPEDSIGNATURES)

Return an array of bools to mark solutions in `res` which fall into `classes` but not `not_classes`
"""
function _get_mask(res, classes, not_classes=[])
    classes == "all" && return fill(trues(length(res.solutions[1])), size(res.solutions))
    bools = vcat([res.classes[c] for c in _vectorise(classes)], [map(.!, res.classes[c]) for c in _vectorise(not_classes)])
    map(.*, bools...)
end


"""
$(TYPEDSIGNATURES)

Go over a solution and an equally-sized array (a "mask") of booleans.
true  -> solution unchanged
false -> changed to NaN (omitted from plotting)
"""
function _apply_mask(solns::Array{Vector{ComplexF64}},  booleans)
    factors = replace.(booleans, 0 => NaN)
    map(.*, solns, factors)
end


# convert x to Float64, raising a warning if complex
function _realify(x)
    !is_real(x) && !isnan(x) ? (@warn "Values with non-negligible complex parts have been projected on the real axis!", x) : nothing
    real(x)
end

_realify(v::Vector) = [_realify.(getindex.(v, k)) for k in 1:length(v[1])]
_realify(a::Array) = _realify.(a)

_vectorise(s::Vector) = s
_vectorise(s) = [s]


# return true if p already has a label for branch index idx
_is_labeled(p::Plots.Plot, idx::Int64) = in(string(idx), [sub[:label] for sub in p.series_list])


function plot1D(res::Result; x::String="default", y::String, class="default", not_class=[], add=false, kwargs...)

    if class == "default"
        if not_class == [] # plot stable full, unstable dashed
            p = plot1D(res; x=x, y=y, class=["physical", "stable"], add=add, kwargs...)
            plot1D(res; x=x, y=y, class="physical", not_class="stable", add=true, style=:dash, kwargs...)
            return p
        else
            p = plot1D(res; x=x, y=y, not_class=not_class, class="physical", add=add, kwargs...)
            return p
        end
    end

    dim(res) != 1 && error("1D plots of not-1D datasets are usually a bad idea.")
    x = x == "default" ? string(first(keys(res.swept_parameters))) : x
    X, Y = transform_solutions(res, [x,y]) # first transform, then filter
    Y = _apply_mask(Y, _get_mask(res, class, not_class))
    branches = _realify(Y)

    # start a new plot if needed
    p = add ? Plots.plot!() : Plots.plot()

    # colouring is matched to branch index - matched across plots
    for k in findall(x -> !all(isnan.(x)), branches[1:end]) # skip NaN branches but keep indices
        l = _is_labeled(p, k) ? nothing : k
        Plots.plot!(_realify.(getindex.(X, k)), branches[k];  color=k, label=nothing, xlabel=latexify(x), ylabel=latexify(y), kwargs...)
    end

    return p
end

plot1D(res::Result, y::String; kwargs...) = plot1D(res; y=y, kwargs...)


function plot2D(res::Result; z::String, branch::Int64, class="physical", not_class=[], add=false, kwargs...)
    X, Y = values(res.swept_parameters)
    Z = getindex.(_apply_mask(transform_solutions(res, z), _get_mask(res, class, not_class)), branch) # first transform, then filter

    p = add ? Plots.plot!() : Plots.plot() # start a new plot if needed

    ylab, xlab = latexify.(string.(keys(res.swept_parameters)))
    p = plot!(map(_realify, [Y, X, Z])...;
    st=:surface, color=:blue, opacity=0.5, xlabel=xlab, ylabel=ylab, zlabel=latexify(z), colorbar=false, kwargs...)
end

function plot2D_cut(res::Result; y::String, cut::Pair, class="default", not_class=[], add=false, kwargs...)

    if class == "default"
        if not_class == [] # plot stable full, unstable dashed
            p = plot2D_cut(res; y=y, cut=cut, class=["physical", "stable"], add=add, kwargs...)
            plot2D_cut(res; y=y, cut=cut, class="physical", not_class="stable", add=true, style=:dash, kwargs...)
            return p
        else
            p = plot2D_cut(res; y=y, cut=cut, not_class=not_class, class="physical", add=add, kwargs...)
            return p
        end
    end

    # the swept params are ranges and thus a sorted search can be performed
    cut_par, cut_value = cut
    cut_par_index = searchsortedfirst(res.swept_parameters[cut_par], cut_value)

    # compare strings beacuse type Num cannot be compared
    swept_pars = res.swept_parameters.keys
    x_index = findfirst(sym -> string(sym)!=string(cut_par), swept_pars)
    isnothing(x_index) && error("The variable $cut_par was not swept over.")
    x = swept_pars[x_index]

    X = res.swept_parameters[x]
    Y =_apply_mask(transform_solutions(res, y), _get_mask(res, class, not_class)) # first transform, then filter
    branches = _realify(x_index==1 ? Y[:, cut_par_index] : Y[cut_par_index, :])

    # start a new plot if needed
    p = add ? Plots.plot!() : Plots.plot()

    # colouring is matched to branch index - matched across plots
    for k in findall(branch -> !all(isnan.(branch)), branches[1:end]) # skip NaN branches but keep indices
        l = _is_labeled(p, k) ? nothing : k
        Plots.plot!(X, branches[k]; color=k, label=l, xlabel=latexify(string(x)), ylabel=latexify(y), kwargs...)
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


plot_phase_diagram(res::Result, class::String; kwargs...) = plot_phase_diagram(res; class=class, kwargs...)


function plot_phase_diagram_2D(res::Result; class="physical", not_class=[], kwargs...)
    X, Y = values(res.swept_parameters)
    Z = sum.(_get_mask(res, class, not_class))

    xlab, ylab = latexify.(string.(keys(res.swept_parameters)))

    # cannot set heatmap ticks (Plots issue #3560)
    heatmap(X, Y, transpose(Z); xlabel=xlab, ylabel=ylab, color=:viridis, kwargs...)
end


function plot_phase_diagram_1D(res::Result; class="physical", not_class=[], kwargs...)
    X = values(res.swept_parameters)
    Y = sum.(_get_mask(res, class, not_class))
    plot(X..., Y; xlabel=latexify(string(keys(res.swept_parameters)...)), ylabel="#", legend=false, yticks=1:maximum(Y), kwargs...)
end

###
# TRANSFORMATIONS TO THE LAB frame
###

function to_lab_frame(soln, res, times)

    timetrace = zeros(length(times))

    for var in res.problem.eom.variables
        val = real(substitute_all(_remove_brackets(var), soln))
        ω = real(substitute_all(var.ω, soln))
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
        ω = real(substitute_all(var.ω, soln))
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







