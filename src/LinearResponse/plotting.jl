function get_jacobian_response(
    res::Result, nat_var::Num, Ω_range, branch::Int; show_progress=true
)
    stable = classify_branch(res, branch, "stable") # boolean array
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")

    spectra = [JacobianSpectrum(res; branch=branch, index=i) for i in findall(stable)]
    C = Array{Float64,2}(undef, length(Ω_range), length(spectra))

    if show_progress
        bar = Progress(
            length(CartesianIndices(C));
            dt=1,
            desc="Diagonalizing the Jacobian for each solution ... ",
            barlen=50,
        )
    end
    # evaluate the Jacobians for the different values of noise frequency Ω
    for ij in CartesianIndices(C)
        C[ij] = abs(evaluate(spectra[ij[2]][nat_var], Ω_range[ij[1]]))
        show_progress ? next!(bar) : nothing
    end
    return C
end
function get_jacobian_response(
    res::Result,
    nat_var::Num,
    Ω_range,
    followed_branches::Vector{Int};
    show_progress=true,
    force=false,
)
    spectra = [
        JacobianSpectrum(res; branch=branch, index=i, force=force) for
        (i, branch) in pairs(followed_branches)
    ]
    C = Array{Float64,2}(undef, length(Ω_range), length(spectra))

    if show_progress
        bar = Progress(
            length(CartesianIndices(C));
            dt=1,
            desc="Diagonalizing the Jacobian for each solution ... ",
            barlen=50,
        )
    end
    # evaluate the Jacobians for the different values of noise frequency Ω
    for ij in CartesianIndices(C)
        C[ij] = abs(evaluate(spectra[ij[2]][nat_var], Ω_range[ij[1]]))
        show_progress ? next!(bar) : nothing
    end
    return C
end

function get_linear_response(
    res::Result, nat_var::Num, Ω_range, branch::Int; order, show_progress=true
)
    stable = classify_branch(res, branch, "stable") # boolean array
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")

    response = ResponseMatrix(res) # the symbolic response matrix
    C = Array{Float64,2}(undef, length(Ω_range), sum(stable))

    # note: this could be optimized by not grabbing the entire huge dictionary every time
    if show_progress
        bar = Progress(
            length(C);
            dt=1,
            desc="Solving the linear response ODE for each solution and input frequency ... ",
            barlen=50,
        )
    end
    for j in findall(stable)

        # get response for each individual point
        s = get_single_solution(res; branch=branch, index=j)
        for i in 1:(size(C)[1])
            C[i, j] = get_response(response, s, Ω_range[i])
        end
        show_progress ? next!(bar) : nothing
    end
    return C
end

function get_rotframe_jacobian_response(
    res::Result, Ω_range, branch::Int; show_progress=true, damping_mod::Float64
)
    stable = classify_branch(res, branch, "stable")
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")
    stableidx = findall(stable)
    C = zeros(length(Ω_range), sum(stable))

    if show_progress
        bar = Progress(
            length(C);
            dt=1,
            desc="Solving the linear response ODE for each solution and input frequency ... ",
            barlen=50,
        )
    end

    for i in 1:sum(stable)
        s = get_single_solution(res; branch=branch, index=stableidx[i])
        jac = res.jacobian(s) #numerical Jacobian
        λs, vs = eigen(jac)
        for j in λs
            for k in 1:(size(C)[1])
                C[k, i] +=
                    1 / sqrt(
                        (imag(j)^2 - Ω_range[k]^2)^2 +
                        Ω_range[k]^2 * damping_mod^2 * real(j)^2,
                    )
            end
        end
        show_progress ? next!(bar) : nothing
    end
    return C
end

"""
    plot_linear_response(res::Result, nat_var::Num; Ω_range, branch::Int, order=1, logscale=false, show_progress=true, kwargs...)

Plot the linear response to white noise of the variable `nat_var` for Result `res` on `branch` for input frequencies `Ω_range`.
Slow-time derivatives up to `order` are kept in the process.

Any kwargs are fed to Plots' gr().

Solutions not belonging to the `physical` class are ignored.
"""
function plot_linear_response(
    res::Result,
    nat_var::Num;
    Ω_range,
    branch::Int,
    order=1,
    logscale=false,
    show_progress=true,
    kwargs...,
)
    length(size(res.solutions)) != 1 &&
        error("1D plots of not-1D datasets are usually a bad idea.")
    stable = classify_branch(res, branch, "stable") # boolean array

    X = Vector{Float64}(collect(values(res.swept_parameters))[1][stable])

    C = if order == 1
        get_jacobian_response(res, nat_var, Ω_range, branch; show_progress=show_progress)
    else
        get_linear_response(
            res, nat_var, Ω_range, branch; order=order, show_progress=show_progress
        )
    end
    C = logscale ? log.(C) : C

    xlabel = latexify(string(first(keys(res.swept_parameters))))
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
function plot_linear_response(
    res::Result,
    nat_var::Num,
    followed_branches::Vector{Int};
    Ω_range,
    logscale=false,
    show_progress=true,
    switch_axis=false,
    force=true,
    kwargs...,
)
    length(size(res.solutions)) != 1 &&
        error("1D plots of not-1D datasets are usually a bad idea.")

    X = Vector{Float64}(collect(first(values(res.swept_parameters))))

    C = get_jacobian_response(res, nat_var, Ω_range, followed_branches; force=force)
    C = logscale ? log.(C) : C

    xlabel = latexify(string(first(keys(res.swept_parameters))))
    ylabel = latexify("Ω")
    if switch_axis
        heatmap(
            Ω_range,
            X,
            C';
            color=:viridis,
            xlabel=ylabel,
            ylabel=xlabel,
            _set_Plots_default...,
            kwargs...,
        )
    else
        heatmap(
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
end

"""
    plot_rotframe_jacobian_response(res::Result; Ω_range, branch::Int, logscale=true, damping_mod::Float64 = 1.0, show_progress=true, kwargs...)

Plot the linear response to white noise in the rotating frame for Result `res` on `branch` for input frequencies `Ω_range`. 'damping_mod' gets multiplied by the real part of the eigenvalues of the Jacobian in order to be able to make peaks with similar frequency seperately identifiable.

Any kwargs are fed to Plots' gr().

Solutions not belonging to the `physical` class are ignored.
"""
function plot_rotframe_jacobian_response(
    res::Result;
    Ω_range,
    branch::Int,
    logscale=true,
    damping_mod::Float64=1.0,
    show_progress=true,
    kwargs...,
)
    length(size(res.solutions)) != 1 &&
        error("1D plots of not-1D datasets are usually a bad idea.")
    stable = classify_branch(res, branch, "stable") # boolean array

    Ω_range = vcat(Ω_range)
    !isempty(findall(x -> x == 0, Ω_range)) &&
        @warn("Probing with Ω=0 may lead to unexpected results")

    X = Vector{Float64}(collect(values(res.swept_parameters))[1][stable])

    C = get_rotframe_jacobian_response(
        res, Ω_range, branch; show_progress=show_progress, damping_mod=damping_mod
    )
    C = logscale ? log.(C) : C

    return heatmap(
        X,
        Ω_range,
        C;
        color=:viridis,
        xlabel=latexify(string(first(keys(res.swept_parameters)))),
        ylabel=latexify("Ω"),
        _set_Plots_default...,
        kwargs...,
    )
end

"""
    plot_eigenvalues(res::Result; branch::Int, class=["physical"], type=:imag, projection=v -> 1, cscheme=:default, kwargs...)

Plot the eigenvalues of the jacobian in the rotating frame for Result `res` on `branch`. Either the real (`type=:real``) or imaginary part (`type=:imag``) can be plotted. The `projection` function ℜᵈ → ℜ is applied to the eigenvectors and defines the color of the eigenvalues. The color scheme can be set to a custom one or to the default one.

Any kwargs are fed to Plots' gr().

Solutions not belonging to the `physical` class are ignored.
"""
function plot_eigenvalues(
    res;
    branch,
    class=["physical"],
    type=:imag,
    projection=v -> 1,
    cscheme=:default,
    kwargs...,
)
    filter = _get_mask(res, class)
    filter_branch = map(x -> getindex(x, branch), replace.(filter, 0 => NaN))

    dim(res) != 1 && error("1D plots of not-1D datasets are usually a bad idea.")
    x = string(first(keys(res.swept_parameters)))
    varied = Vector{Float64}(collect(first(values(res.swept_parameters))))

    eigenvalues = [
        eigvals(res.jacobian(get_single_solution(res; branch=branch, index=i))) for
        i in eachindex(varied)
    ]
    eigenvalues_filtered = map(.*, eigenvalues, filter_branch)

    eigenvectors = [
        eigvecs(res.jacobian(get_single_solution(res; branch=branch, index=i))) for
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
        xlab=latexify(x),
        ylab=ylab,
        zcolor=norm',
        c=myscheme,
        colorbar=false,
        kwargs...,
    )
end
