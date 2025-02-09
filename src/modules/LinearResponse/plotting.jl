"""
$(TYPEDSIGNATURES)

Calculate the Jacobian response spectrum for a given system. Computes the magnitude of the Jacobian response for stable solutions across specified frequency ranges.

# Arguments
- `res::Result`: Result object containing the system's solutions
- `nat_var::Num`: Natural variable to evaluate in the response
- `Ω_range`: Range of frequencies to evaluate
- `branch::Int` or `followed_branches::Vector{Int}`: Branch number(s) to analyze
- `show_progress=true`: Whether to show a progress bar
- `force=false`: Force recalculation of spectrum even if already exists

# Returns
- Array{P,2}: Complex response matrix where rows correspond to frequencies and columns to solutions
"""
function get_jacobian_response(
    res::Result{D,S,P}, nat_var::Num, Ω_range, branch::Int; show_progress=true
) where {D,S,P}
    stable = get_class(res, branch, "stable") # boolean array
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")

    spectra = [JacobianSpectrum(res; branch=branch, index=i) for i in findall(stable)]
    C = Array{P,2}(undef, length(Ω_range), length(spectra))

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
    res::Result{D,S,P},
    nat_var::Num,
    Ω_range,
    followed_branches::Vector{Int};
    show_progress=true,
    force=false,
) where {D,S,P}
    spectra = [
        JacobianSpectrum(res; branch=branch, index=i, force=force) for
        (i, branch) in pairs(followed_branches)
    ]
    C = Array{P,2}(undef, length(Ω_range), length(spectra))

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

"""
$(TYPEDSIGNATURES)

Calculate the linear response of the system for a given branch. Evaluates the linear response by solving the linear response ODE for each stable solution
and input frequency in the given range.

# Arguments
- `res`: Result object containing the system's solutions
- `nat_var::Num`: Natural variable to evaluate in the response
- `Ω_range`: Range of frequencies to evaluate
- `branch::Int`: Branch number to analyze
- `order`: Order of the response to calculate
- `show_progress=true`: Whether to show a progress bar

# Returns
- Array{P,2}: Response matrix where rows correspond to frequencies and columns to stable solutions
"""
function get_linear_response(
    res::Result{D,S,P}, nat_var::Num, Ω_range, branch::Int; order, show_progress=true
) where {D,S,P}
    stable = get_class(res, branch, "stable") # boolean array
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")

    response = ResponseMatrix(res) # the symbolic response matrix
    C = Array{P,2}(undef, length(Ω_range), sum(stable))

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

"""
$(TYPEDSIGNATURES)

Calculate the rotating frame Jacobian response for a given branch. Computes the rotating frame Jacobian response by evaluating eigenvalues of the numerical
Jacobian and calculating the response magnitude for each frequency in the range.

# Arguments
- `res::Result`: Result object containing the system's solutions
- `Ω_range`: Range of frequencies to evaluate
- `branch::Int`: Branch number to analyze
- `show_progress=true`: Whether to show a progress bar
- `damping_mod`: Damping modification parameter

# Returns
- Array{P,2}: Response matrix in the rotating frame

"""
function get_rotframe_jacobian_response(
    res::Result{D,S,P}, Ω_range, branch::Int; show_progress=true, damping_mod
) where {D,S,P}
    stable = get_class(res, branch, "stable")
    !any(stable) && error("Cannot generate a spectrum - no stable solutions!")
    stableidx = findall(stable)
    C = zeros(P, length(Ω_range), sum(stable))

    if show_progress
        bar = Progress(
            length(C);
            dt=1,
            desc="Solving the linear response ODE for each solution and input frequency ... ",
            barlen=50,
        )
    end

    for i in 1:sum(stable)
        s = get_variable_solutions(res; branch=branch, index=stableidx[i])
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
