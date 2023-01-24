using Plots, Latexify, ProgressMeter
export plot_linear_response


function get_jacobian_response(res::Result, nat_var::Num, Ω_range, branch::Int; show_progress=true)
    stable = classify_branch(res, branch, "stable") # boolean array
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")

    spectra = [JacobianSpectrum(res, branch=branch, index = i) for i in findall(stable)]
    C = Array{Float64, 2}(undef,  length(Ω_range), length(spectra))

    if show_progress
        bar = Progress(length(CartesianIndices(C)), 1, "Diagonalizing the Jacobian for each solution ... ", 50)
    end
    # evaluate the Jacobians for the different values of noise frequency Ω
    for ij in CartesianIndices(C)
        C[ij] = abs(evaluate(spectra[ij[2]][nat_var], Ω_range[ij[1]]))
        show_progress ? next!(bar) : nothing
    end
    C
end


function get_linear_response(res::Result, nat_var::Num, Ω_range, branch::Int; order, show_progress=true)

    stable = classify_branch(res, branch, "stable") # boolean array
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")

    response = ResponseMatrix(res) # the symbolic response matrix
    C = Array{Float64, 2}(undef,  length(Ω_range), sum(stable))

    # note: this could be optimized by not grabbing the entire huge dictionary every time
    if show_progress
        bar = Progress(length(C), 1, "Solving the linear response ODE for each solution and input frequency ... ", 50)
    end
    for j in findall(stable)

         # get response for each individual point
        s = get_single_solution(res, branch=branch, index=j)
        for i in 1:(size(C)[1])
            C[i,j] = get_response(response, s, Ω_range[i])
        end
        show_progress ? next!(bar) : nothing
    end
    C
end

function get_rotframe_jacobian_response(res::Result, Ω_range, branch::Int; damping_mod::Float64)
    stable = classify_branch(res, branch, "stable")
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")
    stableidx = findall(stable)
    C = zeros(length(Ω_range), sum(stable));
    print(sum(stable))
    for i in 1:sum(stable)
        s = get_single_solution(res, branch = branch , index = stableidx[i]);
        print(i)
        print("/n")
        jac  = res.jacobian(s) #numerical Jacobian
        λs, vs = eigen(jac)
        for j in λs
            for k in 1:(size(C)[1])
                C[k,i] += 1/sqrt((imag(j)^2-Ω_range[k]^2)^2+Ω_range[k]^2*damping_mod^2*real(j)^2)
            end
        end
        
    end
    C
end


"""
    plot_linear_response(res::Result, nat_var::Num; Ω_range, branch::Int, order=1, logscale=false, show_progress=true, kwargs...)

Plot the linear response to white noise of the variable `nat_var` for Result `res` on `branch` for input frequencies `Ω_range`.
Slow-time derivatives up to `order` are kept in the process.

Any kwargs are fed to Plots' gr().

Solutions not belonging to the `physical` class are ignored.
"""
function plot_linear_response(res::Result, nat_var::Num; Ω_range, branch::Int, order=1, logscale=false, show_progress=true, kwargs...)

length(size(res.solutions)) != 1 && error("1D plots of not-1D datasets are usually a bad idea.")
stable = classify_branch(res, branch, "stable") # boolean array

X = Vector{Float64}(collect(values(res.swept_parameters))[1][stable])

C = order == 1 ? get_jacobian_response(res, nat_var, Ω_range, branch, show_progress=show_progress) : get_linear_response(res, nat_var, Ω_range, branch; order=order, show_progress=show_progress)
C = logscale ? log.(C) : C

heatmap(X, Ω_range,  C; color=:viridis,
    xlabel=latexify(string(first(keys(res.swept_parameters)))), ylabel=latexify("Ω"), HarmonicBalance._set_Plots_default..., kwargs...)
end

function plot_rotframe_jacobian_response(res::Result, Ω_range, branch::Int; logscale=true, damping_mod::Float64 = 1.0, kwargs...)

    length(size(res.solutions)) != 1 && error("1D plots of not-1D datasets are usually a bad idea.")
    stable = classify_branch(res, branch, "stable") # boolean array

    Ω_range = vcat(Ω_range)
    !isempty(findall(x->x==0, Ω_range)) && @warn("Probing with Ω=0 may lead to unexpected results")

    X = Vector{Float64}(collect(values(res.swept_parameters))[1][stable])
    
    C = get_rotframe_jacobian_response(res, Ω_range, branch, damping_mod=damping_mod)
    C = logscale ? log.(C) : C
    
    heatmap(X, Ω_range,  C; color=:viridis,
       xlabel=latexify(string(first(keys(res.swept_parameters)))), ylabel=latexify("Ω"), HarmonicBalance._set_Plots_default..., kwargs...)
end