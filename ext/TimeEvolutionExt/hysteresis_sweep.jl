export plot_1D_solutions_branch


"""Calculate distance between a given state and a stable branch"""
function _closest_branch_index(res::Result,state::Vector{Float64},index::Int64)
    #search only among stable solutions
    physical     = filter_solutions.(
                            filter_solutions.(
                            res.solutions, res.classes["physical"]),res.classes["stable"])

    steadystates = reduce(hcat,physical[index]) #change structure for easier indexing
    distances    = vec(sum(abs2.(steadystates .- state),dims=1))
    return argmin(replace!(distances, NaN=>Inf))
end


"""Return the indexes and values following stable branches along a 1D sweep.
When a no stable solutions are found (e.g. in a bifurcation), the next stable solution is calculated by time evolving the previous solution (quench).

Keyword arguments
- `y`:  Dependent variable expression (parsed into Symbolics.jl) to evaluate the followed solution branches on .
- `sweep`: Direction for the sweeping of solutions. A `right` (`left`) sweep proceeds from the first (last) solution, ordered as the sweeping parameter.
- `tf`: time to reach steady
- `ϵ`: small random perturbation applied to quenched solution, in a bifurcation in order to favour convergence in cases where multiple solutions are identically accessible (e.g. symmetry breaking into two equal amplitude states)
"""
function follow_branch(starting_branch::Int64,res::Result; y="sqrt(u1^2+v1^2)",sweep="right",tf=10000,ϵ=1E-4)
    sweep_directions = ["left", "right"]
    sweep ∈ sweep_directions  || error("Only the following (1D) sweeping directions are allowed:  ", sweep_directions)

    #filter stable solutions
    Y  = transform_solutions(res, y)
    Ys = filter_solutions.(filter_solutions.(Y, res.classes["physical"]),res.classes["stable"]);
    Ys = real.(reduce(hcat,Ys)) #facilitate indexing. Ys is an array (length(sweep))x(#solutions)

    followed_branch    = zeros(Int64,length(Y)) #followed branch indexes
    followed_branch[1] = starting_branch

    if sweep == "left"
        Ys = reverse(Ys,dims=2)
    end

    p1 = collect(keys(res.swept_parameters))[1] #parameter values

    for i in 2:size(Ys)[2]
        s = Ys[followed_branch[i-1],i] #solution amplitude in the current branch and current parameter index
        if !isnan(s) #the solution is not unstable or unphysical
            followed_branch[i] = followed_branch[i-1]
        else #bifurcation found
            if sweep == "right"
                next_index = i
            elseif sweep == "left" #reverse order of iterators
                next_index = length(Y)-i+1
            end

            # create a synthetic starting point out of an unphysical solution: quench and time evolve
            #the actual solution is complex there, i.e. non physical. Take real part for the quench.
            sol_dict = get_single_solution(res, branch=followed_branch[i-1], index= next_index)
            print("bifurcation @ ", p1 ," = ", real(sol_dict[p1])," ")
            sol_dict   = Dict(zip(keys(sol_dict),real.(values(sol_dict)) .+ 0.0im .+ ϵ*rand(length(values(sol_dict)))))
            problem_t  = TimeEvolution.ODEProblem(res.problem.eom, steady_solution=sol_dict, timespan=(0,tf))
            res_t      = TimeEvolution.solve(problem_t,saveat=tf)

            followed_branch[i] = _closest_branch_index(res,res_t.u[end],next_index) #closest branch to final state

            print("switched branch ", followed_branch[i-1] ," -> ", followed_branch[i],"\n")
        end
    end

    if sweep == "left"
        Ys = reverse(Ys,dims=2)
        followed_branch = reverse(followed_branch)
    end

    return followed_branch,Ys
end


"""1D plot with the followed branch highlighted"""
function plot_1D_solutions_branch(starting_branch::Int64,res::Result; x::String, y::String, sweep="right",tf=10000,ϵ=1E-4, xscale="linear",yscale="linear",plot_only=["physical"],marker_classification="stable",filename=nothing,kwargs...)
    plot_1D_solutions(res; x=x, y=y, xscale=xscale,yscale=yscale,
                        plot_only=plot_only,marker_classification=marker_classification,filename=filename,kwargs...)

    followed_branch,Ys = follow_branch(starting_branch,res,y=y,sweep=sweep,tf=tf,ϵ=ϵ)
    if sweep == "left"
        m = "<"
    elseif sweep == "right"
        m = ">"
    end

    Y_followed = [Ys[branch,param_idx] for (param_idx,branch) in enumerate(followed_branch)]

    lines = plot(transform_solutions(res, x),Y_followed,c="k",marker=m; _set_Plots_default..., kwargs...)

    #extract plotted data and return it
    xdata,ydata = [line.get_xdata() for line in lines], [line.get_ydata() for line in lines]
    markers = [line.get_marker() for line in lines]
    save_dict= Dict(string(x) => xdata,string(y)=>ydata,"markers"=>markers)

    !isnothing(filename) ? JLD2.save(_jld2_name(filename), save_dict) : nothing
    return save_dict
end
