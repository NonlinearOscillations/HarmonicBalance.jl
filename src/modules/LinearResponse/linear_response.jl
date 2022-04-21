export get_response
export plot_response


"Get the response matrix corresponding to `res`.
Any substitution rules not specified in `res` can be supplied in `rules`."
function ResponseMatrix(res::Result; rules=Dict())

    # get the symbolic response matrix
    @variables Δ
    M = get_response_matrix(res.problem.eom.natural_equation, Δ, order=2)
    M = substitute_all(M, merge(res.fixed_parameters, rules))
    symbols = cat(res.problem.variables, collect(keys(res.swept_parameters)), dims=1)
    compiled_M = [build_function(el, cat(symbols, [Δ], dims=1)) for el in M]

    ωs = [var.ω for var in res.problem.eom.variables]

    ResponseMatrix(eval.(compiled_M), symbols, ωs)
end


"""Evaluate the response matrix `resp` for the steady state `s` at (lab-frame) frequency `Ω`."""
function evaluate_response_matrix(resp::ResponseMatrix, s::StateDict, Ω)
    values = cat([s[var] for var in resp.symbols], [Ω], dims=1)
    f = resp.matrix
    return [Base.invokelatest(el, values) for el in f]
end


function _evaluate_response_vector(rmat::ResponseMatrix, s::StateDict, Ω)
    m=evaluate_response_matrix(rmat, s, Ω)
    force_pert = cat([[1.0, 1.0*im] for n in 1:size(m)[1]/2]..., dims=1)
    return inv(m) * force_pert
end

"""
$(TYPEDSIGNATURES)

For `rmat` and a solution dictionary `s`,
calculate the total response to a perturbative force at frequency `Ω`.

"""
function get_response(rmat::ResponseMatrix, s::StateDict, Ω)
    resp = 0
    for (i,ω) in enumerate(rmat.ωs)
        this_ω = Float64(substitute_all(ω, s))
        uv1 = _evaluate_response_vector(rmat, s, Ω-this_ω)[2*i-1:2*i]
        uv2 = _evaluate_response_vector(rmat, s, -Ω+this_ω)[2*i-1:2*i]
        resp += sqrt(_plusamp(uv1)^2 + _minusamp(uv2)^2)
    end
    resp
end

# formulas to obtain up- and down- converted frequency components when going from the
# rotating frame into the lab frame
_plusamp(uv) = norm(uv)^2 - 2*(imag(uv[1])*real(uv[2]) - real(uv[1])*imag(uv[2]))
_minusamp(uv) = norm(uv)^2 + 2*(imag(uv[1])*real(uv[2]) - real(uv[1])*imag(uv[2]))


"""
$(TYPEDSIGNATURES)

Plot the linear response of a solution `branch` of `res` for the frequencies `Ω_range`.

This method takes n^2 matrix inversions for n elements of `Ω_range`. It therefore takes longer than `plot_jacobian_spectrum` (which is (O(n)))
but is also valid far from resonance. 

"""
function plot_response(res::Result, Ω_range; branch::Int, logscale=false)
    _set_plotting_settings()
    length(size(res.solutions)) != 1 && error("1D plots of not-1D datasets are usually a bad idea.")
    stability = classify_branch(res, branch, "stable") # boolean array
    !any(stability) && error("Cannot generate a spectrum - no stable solutions!")

    X = Vector{Float64}(collect(values(res.swept_parameters))[1][stability])

    response = ResponseMatrix(res)

    Y = Array{Float64, 2}(undef,  length(Ω_range), length(X))

    # this could be optimized by not grabbing the entire huge dictionary every time

    bar = Progress(length(Y), 1, "Solving the linear response ODE for each solution and input frequency ... ", 50)
    for j in 1:(size(Y)[2])
        next!(bar)
        s = get_single_solution(res, branch=branch, index=j)
        for i in 1:(size(Y)[1])
            Y[i,j] = get_response(response, s, reverse(Ω_range)[i])
        end
    end

    Y = logscale ? log.(Y) : Y
    plt = PyPlot.imshow(Y, aspect="auto", extent=[X[1], X[end], Ω_range[1], Ω_range[end]])
    xlabel(Latexify.latexify(string(first(keys(res.swept_parameters)))),fontsize=20); ylabel("noise angular frequency",fontsize=20);
    return plt
end




