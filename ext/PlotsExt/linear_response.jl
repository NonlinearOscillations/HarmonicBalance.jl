
"""
$(TYPEDSIGNATURES)

Plot the linear response to white noise of the variable `nat_var` for [Result](@ref) `res`
on `branch` identifier.

# Keyword arguments
- `Ω_range`: Range of frequency of the noise probe
- `order`: Order of slow-time derivatives to keep (default: 1)
- `logscale`: Whether to plot response in log scale (default: false)
- `show_progress`: Show progress bar during computation (default: true)
- `kwargs...`: Additional arguments passed to Plots.heatmap

# Returns
A Plots.jl heatmap showing the linear response magnitude across parameter and frequency space.

"""
function HarmonicBalance.plot_linear_response(
    res::Result{D},
    nat_var::Num,
    branch::Int;
    Ω_range,
    order::Int=1,
    logscale::Bool=false,
    show_progress::Bool=true,
    kwargs...,
) where {D}
    D != 1 && error("The results are two dimensional. Consider using the `cut` keyword.")
    stable = HarmonicBalance.get_class(res, branch, "stable") # boolean array

    X = swept_parameters(res)[stable]

    C = if order == 1
        get_jacobian_response(res, nat_var, Ω_range, branch; show_progress)
    else
        get_linear_response(res, nat_var, Ω_range, branch; order=order, show_progress)
    end
    C = logscale ? log.(C) : C

    xlabel = get_labels(res)
    ylabel = latexify("Ω")
    return heatmap(
        X,
        Ω_range,
        C;
        color=:viridis,
        xlabel=xlabel,
        ylabel=ylabel,
        _set_Plots_default...,
        kwargs...,
    )
end

"""
$(TYPEDSIGNATURES)

Plot the linear response to white noise of the variable `nat_var` for [Result](@ref) `res`
on the `followed_branches` identifiers with the size of `Ω_range`.

# Keyword arguments
- `Ω_range`: Range of frequency of the noise probe
- `order`: Order of slow-time derivatives to keep (default: 1)
- `logscale`: Whether to plot response in log scale (default: false)
- `show_progress`: Show progress bar during computation (default: true)
- `kwargs...`: Additional arguments passed to Plots.heatmap

# Returns
A Plots.jl heatmap showing the linear response magnitude across parameter and frequency space.

"""
function HarmonicBalance.plot_linear_response(
    res::Result,
    nat_var::Num,
    followed_branches::Vector{Int};
    Ω_range,
    logscale::Bool=false,
    show_progress::Bool=true,
    switch_axis::Bool=false,
    force::Bool=true,
    kwargs...,
)
    length(size(res.solutions)) != 1 &&
        error("The results are two dimensional. Consider using the `cut` keyword.")

    X = swept_parameters(res)

    C = get_jacobian_response(res, nat_var, Ω_range, followed_branches; force=force)
    C = logscale ? log.(C) : C

    xlabel = get_labels(res)
    ylabel = latexify("Ω")
    if switch_axis
        heatmap(
            Ω_range, X, C'; color=:viridis, ylabel, xlabel, _set_Plots_default..., kwargs...
        )
    else
        heatmap(
            X, Ω_range, C; color=:viridis, xlabel, ylabel, _set_Plots_default..., kwargs...
        )
    end
end

"""
$(TYPEDSIGNATURES)

Plot the linear response to white noise in the rotating frame defined the harmonic ansatz
for [Result](@ref) `res` on `branch` identifier.

# Keyword arguments
- `Ω_range`: Range of frequencies to analyze
- `logscale`: Whether to plot response in log scale (default: true)
- `damping_mod`: Multiplier for the real part of Jacobian eigenvalues (default: 1.0)
- `show_progress`: Show progress bar during computation (default: true)
- `kwargs...`: Additional arguments passed to Plots.heatmap

# Returns
A Plots.jl heatmap showing the response magnitude in the rotating frame.

# Notes
- Setting `damping_mod` < 1 can help distinguish between peaks with similar frequencies
- Solutions not belonging to the `physical` class are ignored
"""
function HarmonicBalance.plot_rotframe_jacobian_response(
    res::Result{D,S,P},
    branch::Int;
    Ω_range,
    logscale=true,
    damping_mod=one(P),
    show_progress=true,
    kwargs...,
) where {D,S,P}
    D != 1 && error("The results are two dimensional. Consider using the `cut` keyword.")
    stable = get_class(res, branch, "stable") # boolean array

    Ω_range = vcat(Ω_range)
    !isempty(findall(x -> x ≈ 0, Ω_range)) &&
        @warn("Probing with Ω=0 may lead to unexpected results")

    X = Vector{P}(collect(values(res.swept_parameters))[1][stable])

    C = get_rotframe_jacobian_response(res, Ω_range, branch; show_progress, damping_mod)
    C = logscale ? log.(C) : C

    return heatmap(
        X,
        Ω_range,
        C;
        color=:viridis,
        xlabel=get_labels(res),
        ylabel=latexify("Ω"),
        _set_Plots_default...,
        kwargs...,
    )
end

"""
$(TYPEDSIGNATURES)

Visualize the eigenvalues of the Jacobian in the rotating frame for `branch` identifier in
the [Result](@ref) `res`.

# Keyword arguments
- `class`: Array of solution classes to include (default: ["physical"])
- `type`: Which part of eigenvalues to plot (`:real` or `:imag`, default: `:imag`)
- `projection`: Function mapping eigenvectors to colors (default: v->1)
- `cscheme`: Color scheme for plotting (`:default` or custom scheme)
- `kwargs...`: Additional arguments passed to Plots.scatter

# Returns
A scatter plot of eigenvalues colored by the projection of their eigenvectors.

# Example
```julia
# Plot imaginary parts of eigenvalues
plot_eigenvalues(result, branch=1)

# Plot real parts with custom coloring based on the norm of eigenvectors of the first harmonic
plot_eigenvalues(result, branch=1, type=:real, projection=v->sqrt(v[1]^2+v[2]^2))
```
"""
function HarmonicBalance.plot_eigenvalues(
    res::Result{D,S,P},
    branch::Int;
    class=["physical"],
    type::Symbol=:imag,
    projection::Function=v -> 1,
    cscheme=:default,
    kwargs...,
) where {D,S,P}
    filter = _get_mask(res, class)
    filter_branch = map(x -> getindex(x, branch), replace.(filter, 0 => NaN))

    D != 1 && error("The results are two dimensional. Consider using the `cut` keyword.")
    varied = Vector{P}(swept_parameters(res))

    eigenvalues = map(eachindex(varied)) do i
        jac = res.jacobian(get_variable_solutions(res; branch=branch, index=i))
        if any(isnan, jac)
            throw(
                ErrorException(
                    "The branch contains NaN values. Likely, the branch has non-physical solutions in the parameter sweep",
                ),
            )
        end
        eigvals(jac)
    end
    eigenvalues_filtered = map(.*, eigenvalues, filter_branch)

    eigenvectors = [
        eigvecs(res.jacobian(get_variable_solutions(res; branch=branch, index=i))) for
        i in eachindex(varied)
    ]
    eigvecs_filtered = map(.*, eigenvectors, filter_branch)

    norm = reduce(
        hcat, [[projection(vec) for vec in eachcol(vecs)] for vecs in eigvecs_filtered]
    )

    if type == :imag
        eigval = reduce(hcat, imag.(eigenvalues_filtered))'
        ylab = L"\Im\{\epsilon\}"
    else
        eigval = reduce(hcat, real.(eigenvalues_filtered))'
        ylab = L"\Re\{\epsilon\}"
    end
    xlab = get_labels(res)

    if cscheme == :default
        colors = theme_palette(cscheme)
        myscheme = cgrad([RGB(1.0, 1.0, 1.0), colors[branch]])
    else
        myscheme = cscheme
    end

    return scatter(
        varied,
        eigval;
        legend=false,
        ms=2,
        markerstrokewidth=0,
        xlab,
        ylab,
        zcolor=norm',
        c=myscheme,
        colorbar=false,
        kwargs...,
    )
end
