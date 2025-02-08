
function phase_diagram(res::Result{D}; kwargs...)
    if D == 1
        plot_phase_diagram_1D(res; _set_Plots_default..., kwargs...)
    elseif D == 2
        plot_phase_diagram_2D(res; _set_Plots_default..., kwargs...)
    else
        error("Data dimension ", D, " not supported")
    end
end

function HarmonicBalance.plot_phase_diagram(res::Result, class::String; kwargs...)
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
    X = first(values(res.swept_parameters))
    Y = sum.(_get_mask(res, class, not_class))
    return plot(
        real.(X),
        Y;
        xlabel=latexify(string(keys(res.swept_parameters)...)),
        ylabel="#",
        legend=false,
        yticks=1:maximum(Y),
        kwargs...,
    )
end
