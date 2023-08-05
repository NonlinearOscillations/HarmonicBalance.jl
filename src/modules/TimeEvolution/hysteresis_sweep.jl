export plot_1D_solutions_branch, follow_branch

"""Calculate distance between a given state and a stable branch"""
function _closest_branch_index(res::Result, state::Vector{Float64}, index::Int64)
    #search only among stable solutions
    stable = HarmonicBalance._apply_mask(res.solutions, HarmonicBalance._get_mask(res, ["physical", "stable"], []))

    steadystates = reduce(hcat, stable[index])
    distances    = vec(sum(abs2.(steadystates .- state), dims=1))
    return argmin(replace(distances, NaN => Inf))
end

"""Return the indexes and values following stable branches along a 1D sweep.
When a no stable solutions are found (e.g. in a bifurcation), the next stable solution is calculated by time evolving the previous solution (quench).
Keyword arguments
- `y`:  Dependent variable expression (parsed into Symbolics.jl) to evaluate the followed solution branches on .
- `sweep`: Direction for the sweeping of solutions. A `right` (`left`) sweep proceeds from the first (last) solution, ordered as the sweeping parameter.
- `tf`: time to reach steady
- `ϵ`: small random perturbation applied to quenched solution, in a bifurcation in order to favour convergence in cases where multiple solutions are identically accessible (e.g. symmetry breaking into two equal amplitude states)
"""
function follow_branch(starting_branch::Int64, res::Result; y="u1^2+v1^2", sweep="right", tf=10000, ϵ=1e-4)
    sweep_directions = ["left", "right"]
    sweep ∈ sweep_directions  || error("Only the following (1D) sweeping directions are allowed:  ", sweep_directions)

    # get stable solutions
    Y  = transform_solutions(res, y)
    Ys = HarmonicBalance._apply_mask(Y, HarmonicBalance._get_mask(res, ["physical", "stable"], []))
    Ys = sweep == "left" ? reverse(Ys) : Ys

    followed_branch    = zeros(Int64,length(Y))  # followed branch indexes
    followed_branch[1] = starting_branch

    p1 = first(keys(res.swept_parameters)) # parameter values

    for i in 2:length(Ys)
        s = Ys[i][followed_branch[i-1]] # solution amplitude in the current branch and current parameter index
        if !isnan(s) # the solution is not unstable or unphysical
            followed_branch[i] = followed_branch[i-1]
        else # bifurcation found
            next_index = sweep == "right" ? i : length(Ys)-i+1

            # create a synthetic starting point out of an unphysical solution: quench and time evolve
            # the actual solution is complex there, i.e. non physical. Take real part for the quench.
            sol_dict  = get_single_solution(res, branch=followed_branch[i-1], index=next_index)

            var = res.problem.variables
            var_values_noise = real.(getindex.(Ref(sol_dict), var)) .+ 0.0im .+ ϵ*rand(length(var))
            for (i, v) in enumerate(var)
                sol_dict[v] = var_values_noise[i]
            end

            problem_t = OrdinaryDiffEq.ODEProblem(res.problem.eom, sol_dict, timespan=(0, tf))
            res_t     = OrdinaryDiffEq.solve(problem_t, OrdinaryDiffEq.Tsit5(), saveat=tf)

            followed_branch[i] = _closest_branch_index(res, res_t.u[end], next_index) # closest branch to final state

            @info "bifurcation @ $p1 = $(real(sol_dict[p1])): switched branch $(followed_branch[i-1]) ➡ $(followed_branch[i])"
        end
    end
    if sweep == "left"
        Ys = reverse(Ys)
        followed_branch = reverse(followed_branch)
    end

    return followed_branch, Ys
end


"""1D plot with the followed branch highlighted"""
function plot_1D_solutions_branch(starting_branch::Int64, res::Result;
    x::String, y::String, sweep="right", tf=10000, ϵ=1e-4, class="default", not_class=[], kwargs...)
    p = plot(res; x=x, y=y, class=class, not_class=not_class, kwargs...)

    followed_branch, Ys = follow_branch(starting_branch, res, y=y, sweep=sweep, tf=tf, ϵ=ϵ)
    Y_followed = [Ys[param_idx][branch] for (param_idx,branch) in enumerate(followed_branch)]
    X = real.(res.swept_parameters[HarmonicBalance._parse_expression(x)])

    Plots.plot!(p, X, real.(Y_followed); linestyle=:dash, c=:gray, label = sweep*" sweep", HarmonicBalance._set_Plots_default..., kwargs...)
    return p
end
