using Plots, Latexify, ProgressMeter
export plot_linear_response


function get_jacobian_response(res::Result, nat_var::Num, Ω_range, branch::Int)
    stable = classify_branch(res, branch, "stable") # boolean array
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")

    spectra = [JacobianSpectrum(res, branch=branch, index = i) for i in findall(stable)]
    C = Array{Float64, 2}(undef,  length(Ω_range), length(spectra))

    bar = Progress(length(CartesianIndices(C)), 1, "Diagonalizing the Jacobian for each solution ... ", 50)
    # evaluate the Jacobians for the different values of noise frequency Ω
    for ij in CartesianIndices(C)
        C[ij] = abs(evaluate(spectra[ij[2]][nat_var], Ω_range[ij[1]]))
        next!(bar)
    end
    C
end


function get_linear_response(res::Result, nat_var::Num, Ω_range, branch::Int; order)

    stable = classify_branch(res, branch, "stable") # boolean array
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")

    response = ResponseMatrix(res) # the symbolic response matrix
    C = Array{Float64, 2}(undef,  length(Ω_range), sum(stable))

    # note: this could be optimized by not grabbing the entire huge dictionary every time
    bar = Progress(length(C), 1, "Solving the linear response ODE for each solution and input frequency ... ", 50)
    for j in findall(stable)

         # get response for each individual point
        s = get_single_solution(res, branch=branch, index=j)
        for i in 1:(size(C)[1])
            C[i,j] = get_response(response, s, Ω_range[i])
        end
        next!(bar)
    end
    C
end


"""
    plot_linear_response(res::Result, nat_var::Num; Ω_range, branch::Int, order=1, logscale=false, kwargs...)

Plot the linear response to white noise of the variable `nat_var` for Result `res` on `branch` for input frequencies `Ω_range`. 
Slow-time derivatives up to `order` are kept in the process.

Any kwargs are fed to Plots' gr().

Solutions not belonging to the `physical` class are ignored. 
"""
function plot_linear_response(res::Result, nat_var::Num; Ω_range, branch::Int, order=1, logscale=false, kwargs...)

length(size(res.solutions)) != 1 && error("1D plots of not-1D datasets are usually a bad idea.")
stable = classify_branch(res, branch, "stable") # boolean array

X = Vector{Float64}(collect(values(res.swept_parameters))[1][stable])

C = order == 1 ? get_jacobian_response(res, nat_var, Ω_range, branch) : get_linear_response(res, nat_var, Ω_range, branch; order=order)
C = logscale ? log.(C) : C

heatmap(X, Ω_range,  C; color=:viridis, 
    xlabel=latexify(string(first(keys(res.swept_parameters)))), ylabel=latexify("Ω"), HarmonicBalance._set_Plots_default..., kwargs...)
end




