using PyPlot
using PyCall
using Latexify
using JLD2
export plot_1D_solutions
export plot_2D_phase_diagram
export plot_2D_phase_diagram_interactive
export transform_solutions
export _set_plotting_settings


"Set global plotting settings"
function _set_plotting_settings()
    plt.style.use("default") #reset settings
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams") 
    rcParams["text.usetex"] = true
    rcParams["font.family"]       = "sans-serif"
    rcParams["font.sans-serif"]   = ["Helvetica"]
    rcParams["legend.frameon"]    = false
    rcParams["legend.fontsize"]   = "x-large"
    rcParams["patch.linewidth"]   =  0.3
    rcParams["lines.linewidth"]   = 2
    rcParams["lines.markersize"]  = 7
    rcParams["grid.linestyle"]    =  "-"
    rcParams["grid.linewidth"]    =  0.75
    rcParams["axes.axisbelow"]    =  true
    rcParams["axes.linewidth"]    =  1.25
    rcParams["axes.labelsize"]    = "large"
    rcParams["xtick.labelsize"]   = "xx-large"
    rcParams["ytick.labelsize"]   = "xx-large"
    rcParams["xtick.direction"]   = "in"
    rcParams["ytick.direction"]   = "in"
    rcParams["xtick.major.size"]  = 6
    rcParams["ytick.major.size"]  = 6
    rcParams["xtick.minor.size"]  = 3
    rcParams["ytick.minor.size"]  = 3
    rcParams["xtick.major.width"] = 1
    rcParams["ytick.major.width"] = 1
    rcParams["xtick.minor.width"] = .5
    rcParams["ytick.minor.width"] = .5
    rcParams["xtick.major.pad"]   = 7
    rcParams["ytick.major.pad"]   = 7
    rcParams["figure.autolayout"] = true
end



# the total number of solutions
total_elements(array) = prod(size(array))

add_dim(x::Array) = reshape(x, (1,size(x)...)) #quickfix for the case where axx is not a matrix. It would require an extra singleton dimension

remove_singleton!(arr) = dropdims(arr, dims = tuple(findall(size(arr) .== 1)...)) #remove singleton dimensions of an array

"""
$(TYPEDSIGNATURES)

Goes over a solution and an equally-sized array (a "mask") of booleans. 
true  -> solution unchanged
false -> changed to NaN (omitted from plotting)
"""
function filter_solutions(solution::Vector,  booleans)
    total_elements(solution) == total_elements(booleans) || error("attempt to filter a solution using a wrongly-sized boolean array")
    rules = Dict(1 => 1., 0 => NaN)
    factors = [rules[pt] for pt in booleans]
    return solution .* factors
end


"""
$(TYPEDSIGNATURES)

Takes a `Result` object and a string `f` representing a Symbolics.jl expression.
Returns an array with the values of `f` evaluated for the respective solutions.
Additional substitution rules can be specified in `rules` in the format `("a" => val)` or `(a => val)`
"""
function transform_solutions(res::Result, f::String; rules=Dict())
    # a string is used as input - a macro would not "see" the user's namespace while the user's namespace does not "see" the variables
    transformed = [Vector{ComplexF64}(undef, length(res.solutions[1])) for k in res.solutions] # preallocate
        
    # define variables in rules in this namespace
    new_keys = declare_variable.(string.(keys(Dict(rules)))) 
    expr = f isa String ? Num(eval(Meta.parse(f))) : f

    fixed_subs = merge(res.fixed_parameters, Dict(zip(new_keys, values(Dict(rules)))))
    expr = substitute_all(expr, Dict(fixed_subs))

    vars = res.problem.variables
    all_symbols = cat(vars, collect(keys(res.swept_parameters)), dims=1)
    comp_func = build_function(expr, all_symbols) 
    f = eval(comp_func)

    # preallocate an array for the numerical values, rewrite parts of it
    # when looping through the solutions
    vals = Vector{ComplexF64}(undef, length(all_symbols))
    n_vars = length(vars)
    n_pars = length(all_symbols) - n_vars

    for idx in CartesianIndices(res.solutions)
        params_values = res.swept_parameters[Tuple(idx)...]
        vals[end-n_pars+1:end] .= params_values # param values are common to all branches
        for (branch,soln) in enumerate(res.solutions[idx])
            vals[1:n_vars] .= soln
            transformed[idx][branch] = Base.invokelatest(f, vals)
        end
    end
    return transformed
end


"Construct a matrix of subplots from a linear subplot array, even if input length and # of rows or columns is not commensurable"
function resize_axes!(f,axs,nrows,ncols)
    gspec = pyimport("matplotlib.gridspec")
    gs = gspec.GridSpec(nrows,ncols,hspace=0.6,wspace=0.6)

    if nrows >= ncols
        for i in 1:nrows
            for j in 0:ncols
                k = i + j*ncols
                if k <= length(axs)
                    axs[k].set_position(gs[k].get_position(f))
                end
            end
        end
    else
        for i in 1:nrows
            for j in 0:ncols
                k = i + j*nrows
                if k <= length(axs)
                    axs[k].set_position(gs[k].get_position(f))
                end
            end
        end
    end
    return axs
end


"""
    plot_1D_solutions(res::Result; 
                        x::String, y::String, 
                        x_scale=1.0, y_scale=1.0, 
                        marker="o",xscale="linear",yscale="linear"
                        ,plot_only=["physical"],
                        marker_classification="stable",filename=nothing)

Make a 1D plot of a `Result` object.    

Keyword arguments
- `x`, `y`: Expressions to plot on as independent/dependent variables (parsed into Symbolics.jl).
- `x_scale`, `y_scale`: Factors to multiply the shown axis ticks with.
- `marker`: The point marker to use.
- `xscale`, `yscale` = x_scale.
- `plot_only`: a list of strings corresponding to the solution classes of `Result`. Only solutions which belong to the listed classes are plotted.
- `marker_classification`: A class of the solutions (created by `classify_solutions!`) which is distinguished with different markers. Entering an inequality creates a new class "custom_class".
- `filename`: if different from `nothing`, plotted data and parameter values are exported to `./filename.jld2`.

The strings in `marker_classification` allows the user to stablish custom criteria for binary classification of solutions. For instance, if `marker_classification = "ω^15* sqrt(u1^2 + v1^2) < 0.1"`, 
for a system with harmonic variables u1,v1, then solutions are classified as `true` according to that criterion and `false` according to its complement. 
"""
function plot_1D_solutions(res::Result; x::String, y::String, x_scale=1.0, y_scale=1.0, marker="o",xscale="linear",yscale="linear",plot_only=["physical"],marker_classification="stable",filename=nothing)
    _set_plotting_settings()
    length(size(res.solutions)) != 1 && error("1D plots of not-1D datasets are usually a bad idea.")

    X = transform_solutions(res, x)
    Y = transform_solutions(res, y) # first transform, then filter
    
    #retrieve data from specific branch
    relevant_indices = findall(x->x==true, [!all(isnan.(getindex.(Y, branch))) for branch in 1:length(Y[1])])
    X = x_scale .* [p[relevant_indices] for p in X]
    Y = y_scale .* [p[relevant_indices] for p in Y]

    f,ax = subplots(1,1,figsize=(6,5))

    #filtering of the solutions according to keyword arguments
    "stable" in plot_only && "physical" ∉ plot_only && error("Stability is not defined for unphysical solutions!")
    ~(sum([class in keys(res.classes) for class in plot_only])==length(plot_only)) && error("Class not found. Solutions are classified according to ", keys(res.classes))
    ("binary_labels" in plot_only) && print("Binary label is not Boolean and thus was ignored in plotting")
    
    for class in plot_only 
        if  class!="binary_labels" 
            Y = filter_solutions.(Y, res.classes[class])
        end
    end

    sol_type, not_sol_type = _classify_plot_data(res, marker_classification)

    Ys, Yu = filter_solutions.(Y, res.classes[sol_type]), filter_solutions.(Y, [.!el for el in res.classes[sol_type]])
    lines = ax.plot(X, Ys, marker) #Nan solutions are ignored 
    ax.set_prop_cycle(nothing) #reset color cycler state (NaN solutions aren't shown but the color cycler runs)
    append!(lines,ax.plot(X, Yu, "X"))

    leg2 = [plt.Line2D([0], [0], marker=marker, color="w", label=sol_type, markerfacecolor="k", markersize=10),
            plt.Line2D([0], [0], marker="X", color="w", label=not_sol_type,markerfacecolor="k", markersize=10)] 

    if !isnothing(filename)
        xdata,ydata = [line.get_xdata() for line in lines], [line.get_ydata() for line in lines]
        markers = [line.get_marker() for line in lines]
        marker_dict = Dict(marker=>sol_type,"X"=>not_sol_type)
        JLD2.save(_jld2_name(filename), Dict(string(x) => xdata,string(y)=>ydata,"marker_dict"=>marker_dict,"markers"=>markers))
    end

    ax.set_xlabel(Latexify.latexify(x),fontsize=24) 
    ax.set_ylabel(Latexify.latexify(y),fontsize=24) 
    leg1 = legend(string.(relevant_indices))
    ax.add_artist(leg1)
    
    ax.legend(handles=leg2,loc="center right") 
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
end


""" For a plotted `marker_classification`, classify the data in `res`. 
If `marker_classification` is an existing class, choose that class. If an inequality, parse it to create a new class."""
function _classify_plot_data(res, marker_classification)
    if marker_classification ∈ keys(res.classes)
        sol_type, not_sol_type = marker_classification, "not " * marker_classification
    elseif ~isnothing(marker_classification) # if an inequality is supplied, create a new class
        try
            classify_solutions!(res, marker_classification, "custom class")
        catch
            error(marker_classification * " is not a valid classification or an inequality!")
        end
        sol_type, not_sol_type = "custom class", "not custom class"
    end
    sol_type, not_sol_type
end




"""
    plot_1D_jacobian_eigenvalues(res::Result; x::String, physical=true, stable=false,marker_re="o",marker_im="X", filename=nothing)

Make a 1D plot of the Jacobian eigenvalues for each of the solutions in a `Result` object.

Keyword arguments

- `x`: The function on the x axis (a string parsed into Symbolics.jl).
- `physical`, `stable`: Booleans specifying whether unphysical and/or unstable solutions are shown.
- `marker_re`, `marker_im`: The markers to use for the Re and Im parts of the eigenvalues.
- `filename`: if different from `nothing`, plotted data and parameter values are exported to `./filename.jld2`.

"""
function plot_1D_jacobian_eigenvalues(res::Result; x::String, physical=true, stable=false,marker_re="o",marker_im="X", filename=nothing)
    _set_plotting_settings()

    xplot = transform_solutions(res, x) #indepedenent variable to plot
    
    #filtering of the solutions according to keyword arguments
    !physical && stable && error("Stability is not defined for unphysical solutions!")
    
    numbers = [Base.not_int.(isnan.(prod.(s))) for s in res.solutions] # true for solution which are not NaN
    to_evaluate = physical ? (stable ? res.classes["stable"] : res.classes["physical"]) : numbers
    to_evaluate = [to_evaluate[i] .* numbers[i] for (i,_) in enumerate(numbers)] 

    nsolsmax  = length(res.solutions[1])
    f,ax = subplots(1,nsolsmax,figsize=(4*nsolsmax,4))
    if nsolsmax==1
        ax = add_dim([ax])
    end

    lines_re,lines_im =[],[]
    for branch in 1:nsolsmax
        λs = [NaN*ones(ComplexF64, length(res.problem.variables)) for i in 1:length(res.solutions)]
        for index in 1:length(res.solutions)
            # calculate if in to_evaluate, otherwise leave the NaN
            if to_evaluate[index][branch]
                J = res.jacobian(get_single_solution(res, branch=branch, index=index))
                λs[index] = eigvals(J)
            end
        end

        X_plot = real.([xplot[idx][branch] for idx in 1:length(res.solutions)])
        append!(lines_re,ax[branch].plot(X_plot,real.(λs),string(marker_re,"r")))
        append!(lines_im,ax[branch].plot(X_plot,imag.(λs),string(marker_im,"g")))  #
        ax[branch].set_title(string("solution ", branch),fontsize=12,loc="left"); 
        ax[branch].set_xlabel(Latexify.latexify(string(x)),fontsize=24)
        ax[branch].axhline(0,ls=":",lw=4) #reference value for unstability test 
    end 

    if !isnothing(filename)
        save_dict =  Dict(zip(["re","im"],[Dict("data"=>Dict(),"marker"=>[]) for d in ["re","im"]]))
        for (axis,lines,marker) in zip(["re","im"],[lines_re,lines_im],[marker_re,marker_im])
            xdata = [line.get_xdata() for line in lines]
            ydata = [line.get_ydata() for line in lines]
            save_dict[axis]["data"] = Dict(string(x) => xdata,string(axis," part(eig)")=>ydata)
            save_dict[axis]["marker"] = marker
        end
        
        JLD2.save(_jld2_name(filename), save_dict)
    end

    legend_elements = [plt.Line2D([0], [0], marker=marker_re, color="w", label=L"\Re({\mathrm{eig}(J)})",
                        markerfacecolor="r", markersize=7),
                        plt.Line2D([0], [0], marker=marker_im, color="w", label=L"\Im({\mathrm{eig}(J)})",
                        markerfacecolor="g", markersize=7)]    
   ax[end].legend(handles=legend_elements,loc="best",fontsize=15) 
end

"""
    plot_2D_solutions(res::Result,ax=nothing; filename=nothing)

Make a 2D plot of each of solutions vs swept parameters for a `Result` object, obtained from a `get_steady_states` applied to a 2D parameter grid.`.

Keyword arguments

- `ax`: axis object from `PyCall.PyObject` setting the coordinate system where data will be plotted. If not given, it is created automatically.
- `filename`: if different from `nothing`, plotted data and parameter values are exported to `./filename.jld2`.
"""
function plot_2D_solutions(res::Result; ax=nothing, filename=nothing)
    _set_plotting_settings()
    nvar  = length(res.solutions[1,1][1]) #number of variables
    nsols = length(res.solutions[1,1]) #maximum number of solutions
    x,y = collect(values(res.swept_parameters))
    extent=[x[1],x[end],y[1],y[end]]

    #transform the solution structs into tensors that are easier to handle by imshow
    physical_solutions = real.(filter_solutions.(res.solutions,res.classes["physical"]))
    physical_sols = reshape(reduce(hcat,[reduce(hcat,sol) for sol in physical_solutions]),nvar,nsols,length(x),length(y))
    
    var_names = _get_var_name_labels(res) #variable names for plot labels

    input_ax = ax
    if isnothing(input_ax) #create figure if axes are not provided, otherwise accept input
        f,ax = subplots(nsols,nvar,figsize=(4*nvar,4*nsols))
    end

    save_dict = Dict([string("panel (",m,",",l,")")=> Dict() for m in 1:nvar for l in 1:nsols])
    px,py = string.(collect(keys(res.swept_parameters)))  #swept parameter strings
    for m in 1:nvar
        for l in 1:nsols 
            a = ax[l,m].imshow(physical_sols[m,l,:,end:-1:1]',extent=extent,aspect="auto")
            ax[l,m].set_xlabel(Latexify.latexify(px),fontsize=24); 
            ax[l,m].set_ylabel(Latexify.latexify(py),fontsize=24); 
            if !isnothing(filename)
                save_dict[string("panel (",m,",",l,")")]= Dict("variable"=>var_names[m],"solution #"=>l,"data"=>a.get_array(),
                string("(",px,"_min ",px,"_max ",py,"_min ",py,"_max)")=>extent)
            end
        end
        ax[1,m].set_title(Latexify.latexify(var_names[m]),fontsize=20)
    end
    f.tight_layout()

    if !isnothing(filename)
        JLD2.save(_jld2_name(filename), save_dict)
    end
end


"""
    plot_2D_phase_diagram(res::Result; stable=false,observable="nsols",ax=nothing, filename=nothing)

Make a 2D phase diagram plot of each of solutions vs swept parameters for a `Result` object, obtained from a `get_steady_states` applied to a 2D parameter grid.

Keyword arguments 

- `stable`: whether only stable solutions are depicted
- `observable`: reference observable to represent dynamical phases in the problem. If `observable="nsols"`, number of solutions for each point is shown. 
   If instead `observable="binary"`, the result of classification of bistrings `[is_stable(solution_1),is_stable(solution_2),...]` is presented (see `classify_binaries!(Result)` function).
- `filename`: if different from `nothing`, plotted data and parameter values are exported to `./filename.jld2`.

"""
function plot_2D_phase_diagram(res::Result; stable=false,observable="nsols",ax=nothing, filename=nothing)
    _set_plotting_settings()
    observable_choice = ["nsols", "binary"]
    observable ∈ observable_choice || error("Only the following 2D observables are allowed:  ", observable_choice)

    X,Y = collect(values(res.swept_parameters))
    x,y = collect(keys(res.swept_parameters))
    extent = [X[1], X[end], Y[1], Y[end]]

    if observable=="nsols" #number of solutions plot
        if stable
            obs_2D =  sum.(res.classes["stable"])
            pd_title = "number of stable solutions"
        else
            obs_2D =  sum.(res.classes["physical"])
            pd_title = "total number of solutions"
        end
    elseif observable=="binary"
        obs_2D = res.classes["binary_labels"]
        pd_title = "binary classified solutions"
    end

    input_ax = ax
    if isnothing(input_ax) #create figure if axes are not provided, otherwise accept input
        f,ax = subplots(1,1,figsize=(6,5))
    end
    # rearrange num_of_sols to match imshow output
    im = ax.imshow(obs_2D[:,end:-1:1]', extent=extent, aspect="auto",vmin=minimum(obs_2D),vmax=maximum(obs_2D))
    if isnothing(input_ax) 
        Nmax = maximum(obs_2D)
        colorbar(im,ax=ax,ticks=collect(1:Nmax), boundaries=collect(1:Nmax+1).-0.5)
    end

    px,py =string.([x,y])#swept parameter strings

    ax.set_xlabel(Latexify.latexify(px),fontsize=24)
    ax.set_ylabel(Latexify.latexify(py),fontsize=24)
    ax.set_title(pd_title)

    if !isnothing(filename)
        JLD2.save(_jld2_name(filename), Dict("observable"=>observable,"data"=>im.get_array(),
                                                 string("(",px,"_min ",px,"_max ",py,"_min ",py,"_max)")=>extent))
    end
end


#################################################
###interactive plotting 
#################################################

"Make 2-D coordinate arrays for vectorized evaluations of 2-D scalar/vector fields over 2-D grids"
function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

"Argument corresponding to the closest value of an array to a given point"
_find_nearest(array,value) = argmin(abs.(array.-value)) 

_reshape_for_plot(arr,L,M,N) = reshape(reduce(hcat,reduce(hcat,arr)),L,M,N) #some useful transformation from matrices of vectors to higher dimensional tensors, for plotting

_get_var_name_labels(res::Result) = string.(res.problem.variables) 

"Preprocess solution arrays to be plotted. Extract parameter information for labelling and setting axes of interactive plots"
function _get_interactive_plot_variables(res::Result,cut_type; string_f, marker_classification)
    x,y = collect(keys(res.swept_parameters))

    X,Y = collect(values(res.swept_parameters))

    gx,gy = collect(meshgrid(X, Y)) #artificial grid points
    var_names = _get_var_name_labels(res) #variable names for plot labels

    nvars     = length(res.solutions[1][1])
    nsolsmax  = length(res.solutions[1])

    sol_type, not_sol_type = _classify_plot_data(res, marker_classification)
    
    #prepare solution arrays for plotting
    if cut_type =="transform"
        solutions = [transform_solutions(res, string) for string in string_f]
        physical_sols = [filter_solutions.(solution, res.classes["physical"]) for solution in solutions]

        Ys = [real.(filter_solutions.(physical_sol, res.classes[sol_type]))  for physical_sol in physical_sols]
        Yu = [real.(filter_solutions.(physical_sol, [.!el for el in res.classes[sol_type]])) for physical_sol in physical_sols]
    else
        physical_sols = filter_solutions.(res.solutions, res.classes["physical"])
        Ys = real.(filter_solutions.(physical_sols, res.classes[sol_type])) 
        Yu = real.(filter_solutions.(physical_sols, [.!el for el in res.classes[sol_type]]))
    end

    return nvars,nsolsmax,Ys,Yu,x,y,X,Y,gx,gy,var_names,sol_type,not_sol_type
end


"Set up axes and plot an invisible grid of points to be used as reference for subsequent hovering labels"
function _get_interactive_plot_axes(x,y,gx,gy,var_names,cut_dim,cut_type,nvars,nsolsmax,sol_type,not_sol_type; string_f)
    if  cut_type=="solutions"    
        N_panels = nvars + 1
        lab = [sol_type,not_sol_type]
    elseif cut_type=="jacobian_eigenvalues"
        N_panels = nsolsmax + 1
        lab = ["real","imag"]
    elseif cut_type == "transform"
        N_panels = length(string_f) + 1
        lab =   [sol_type,not_sol_type]  
    end    

    f,ax = plt.subplots(1,N_panels,figsize=(5*N_panels,5))

    sc = ax[1].scatter(gx',gy', s=1,alpha=0.) #create an invisible grid to attribute hovering labels, the transposes are introduce to match imshow ordering
    ax[1].set_xlabel(Latexify.latexify(x),fontsize=24); ax[1].set_ylabel(Latexify.latexify(y),fontsize=24)

    #hovering label preparation
    annot = ax[1].annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                    bbox=Dict(("boxstyle"=>"round"), ("fc"=>"w")),
                    arrowprops=Dict("arrowstyle"=>"->"))
    annot.set_visible(false)

    for l in 1:N_panels-1 #preparation of plot labels, no need to update these onclick
        if cut_dim=="1" #select direction along which solution cut_dim will be drawn
            ax[l+1].set_xlabel(Latexify.latexify(string(x)),fontsize=24); 
        elseif cut_dim=="2"
            ax[l+1].set_xlabel(Latexify.latexify(string(y)),fontsize=24); 
        end

        if cut_type=="solutions"
            ax[l+1].set_ylabel(Latexify.latexify(var_names[l]),fontsize=24);
        elseif cut_type=="jacobian_eigenvalues" 
            ax[l+1].set_title(string("solution ", l),fontsize=12,loc="left"); 
        elseif cut_type=="transform"
            ax[l+1].set_title(string_f[l])
        end    
    end
    
    return sc,ax,f,annot,lab,im
end


"Plotter of Jacobian eigenvalues along a cut of a 2D solutions. TODO: condense this with  plot_1D_jacobian_eigenvalues"
function _plot_2D_solutions_jacobian_cut(res::Result,filtered_sol,parameter_cut,Z::Vector{Float64},ax::PyCall.PyObject)
    
     # pre-evaluate the Jacobian with fixed parameters
    fixed_pair     = parameter_cut[2]
    fixed_params = Dict(k=>parse(ComplexF64,string(v))  for (k,v) in pairs(merge(res.fixed_parameters,Dict(fixed_pair))))
    fixed_params = Num_to_Variable(fixed_params)
    
    Jac = HomotopyContinuation.evaluate(res.problem.jacobian,collect(keys(fixed_params))=>collect(values(fixed_params))) #maybe make this work with HarmonicBalance.Jacobian?
    
    swept_p      = Num_to_Variable(parameter_cut[1][1]) #swept parameter symbolic variable
    swept_values = parameter_cut[1][2] #swept parameter values

     # evaluate the numerical jacobian at the solutions. NaN elements in solutions need to be removed. These are introduced by filter_solutionss.
    Js = [HomotopyContinuation.evaluate(Jac, 
            res.problem.system.variables => filtered_sol[idz,:],
            swept_p=> swept_values[idz]) for idz in 1:length(Z) if any(isnan.(filtered_sol[idz,:]).==false)]
   
    evals = [eigvals(J) for J in Js]
    Z_plot = [Z[idz] for idz in 1:length(Z) if any(isnan.(filtered_sol[idz,:]).==false)]
    ax.plot(Z_plot,real.(evals),"-")
    ax.plot(Z_plot,imag.(evals),"--")
end



"Take cut of filtered solutions once the user clicks on a 2D plot"
function _prepare_solution_cuts(ax::PyCall.PyObject,res::Result,Ys,Yu,cut_dim,cut_type,idx,idy,nvars,nsolsmax,X,Y)
    p1,p2 = collect(keys(res.swept_parameters))
    if cut_dim=="1"  #select direction along which solution cut_dim will be drawn
        ax.axhline(Y[idy],ls="--",c="w")   
        if cut_type=="transform"
            solution_cut_s   = reduce(hcat,Ys[:,idy])
            solution_cut_u   = reduce(hcat,Yu[:,idy])
        else
            solution_cut_s   = _reshape_for_plot(Ys[:,idy],nvars,nsolsmax,length(X)) 
            solution_cut_u   = _reshape_for_plot(Yu[:,idy],nvars,nsolsmax,length(X)) 
        end
        parameter_val  = [(p1=>res.swept_parameters[p1]), (p2=>res.swept_parameters[idy][2])]
        Z = X #cut parameter
    elseif cut_dim=="2"  
        ax.axvline(X[idx],ls="--",c="w")
        if cut_type=="transform"
            solution_cut_s   = reduce(hcat,Ys[idx,:])
            solution_cut_u   = reduce(hcat,Yu[idx,:])
        else
            solution_cut_s  =  _reshape_for_plot(Ys[idx,:],nvars,nsolsmax,length(Y)) 
            solution_cut_u  =  _reshape_for_plot(Yu[idx,:],nvars,nsolsmax,length(Y))   
        end
        parameter_val  = [(p2=>res.swept_parameters[p2]), (p1=>res.swept_parameters[idx][1])]
        Z = Y #cut parameter  
    end    
    return Z,solution_cut_s, solution_cut_u,parameter_val      
end

"""
    plot_2D_phase_diagram_interactive(res::Result;
     observable="nsols",
      stable=false,nrows=2,
      ncols=2,cut_dim="1",
      cut_type="solutions",
      string_f=nothing,
      marker_classification="stable")

Interactive phase diagram of 2D solutions stored in `Result`. 
This includes a clickable version of a `plot_2D_phase_diagram` for a given `observable` and extra panels containing 1D cuts of solutions, functions of solutions, or Jacobian eigenvalues.
      

Keyword arguments 

- `observable`: reference observable to represent dynamical phases in the problem. If `observable="nsols"`, number of solutions for each point is shown. 
   If instead `observable="binary"``, the result of classification of bistrings `[is_stable(solution_1),is_stable(solution_2),...]` is presented (see `classify_binaries!(Result)` function).
- `ncols`, `nrows`: number of rows and columns of the plot window.
- `cut_dim`: dimension along which 1D quantities will be calculated. `cut_dim="1"` (`cut_dim="2"`) takes a cut along the horizontal (vertical) parameter dimension of the 2D plot.
- `cut_type`: quantity to be represented along the 1D cut. If `cut_type=solutions`, steady state variables are shown with a panel per `Problem` variable.
    Else if `cut_type=jacobian eigenvalues`, Re and Im parts of complex Jacobian eigenvalues for each solution are shown  with a panel per solution. 
    If instead `cut_type=transform`, functions of the solution variables passed to `string_f` (see below) are displayed.
- `string_f`: list of strings for transformed observables to be plotted in 1D when `cut_type=transform`, e.g. `string_f=["sqrt(u1^2 + v1^2)","sqrt(u2^2 + v2^2)"]` for `Problem` variables `u1,u2,v1,v2`.
- `marker_classification`: A class of the solutions (created by `classify_solutions!`) which is distinguished with different markers. Entering an inequality creates a new class "custom_class".

"""
function plot_2D_phase_diagram_interactive(res::Result; observable="nsols", stable=false,nrows=2,ncols=2,cut_dim="1",cut_type="solutions",string_f=nothing,marker_classification="stable")
    pygui(true) #opens a separate window
    _set_plotting_settings()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams") 
    rcParams["text.usetex"] = false #cannot be set globally as it seems to conflict with pygui(true)

    #label each location in parameter space by string stating which solutions are stable/unstable. Consider only physical solutions
    names=vec(join.(replace!.([[join.(string(el)) for el in bit_string[phys_string]] for (bit_string,phys_string) in zip(res.classes["stable"],res.classes["physical"])], ("true"=>"S"),("false"=>"U"))));
    
    cut_types = ["solutions", "jacobian_eigenvalues","transform"]
    cut_type ∈ cut_types || error("Only the following types of 1D cuts are allowed:  ", cut_types)
    
    nvars,nsolsmax,Ys,Yu,x,y,X,Y,gx,gy,var_names,sol_type,not_sol_type = _get_interactive_plot_variables(res,cut_type,string_f=string_f,marker_classification=marker_classification)
    sc,ax,f,annot,lab,im = _get_interactive_plot_axes(x,y,gx,gy,var_names,cut_dim,cut_type,nvars,nsolsmax,sol_type,not_sol_type,string_f=string_f)
    
    length(vec(ax)) <= nrows*ncols || error("insufficient # of panels requested, please increase nrows or ncols") #sanity check before any plot is made
   
    im = plot_2D_phase_diagram(res; stable=stable,observable=observable,ax=ax[1])

    "Update annotations when cursor moves"
    function update_annot(ind) 
        pos = sc.get_offsets()[ind["ind"][1],:]
        annot.xy = pos
        text = names[ind["ind"][1]]

        annot.set_text(text)
        annot.get_bbox_patch().set_facecolor("w")
        annot.get_bbox_patch().set_alpha(0.8)
    end

    "Action when mouse hovers over the figure: show string indicating stability/instability of solution"
    function hover(event) 
        vis = annot.get_visible()
        if event[:inaxes] == ax[1]
            cont, ind = sc.contains(event)
            if cont
                update_annot(ind)
                annot.set_visible(true)
                f[:canvas].draw_idle()
            else
                if vis
                    annot.set_visible(false)
                    f[:canvas].draw_idle()
                end
            end
        end
    end
               
    f[:canvas][:mpl_connect]("motion_notify_event",hover)    
     

    "Simple mouse click function to store coordinates and plot the corresponding cuts"
    function onclick(event)           
        ix, iy = event[:xdata], event[:ydata]
        idx = _find_nearest(X, ix)  #find closest point to the mouse click.  
        idy = _find_nearest(Y, iy)     
       
        if length(ax[1][:lines])>0 #clear plot in current axis in case there is any
            ax[1][:lines] = []
        end

        #Z is the parameter along which cut is taken
        if cut_type!="transform"  
            Z,solution_cut_s, solution_cut_u,parameter_val  = _prepare_solution_cuts(ax[1],res,Ys,Yu,cut_dim,cut_type,idx,idy,nvars,nsolsmax,X,Y) #ax[1] is the axis where obs_2D is displayed
        else
            solution_cut_s   = []
            solution_cut_u = []
            for l=1:length(string_f)
                Z,sol_cut, sol_cut_u,parameter_val  = _prepare_solution_cuts(ax[1],res,Ys[l],Yu[l],cut_dim,cut_type,idx,idy,nvars,nsolsmax,X,Y)
                append!(solution_cut_s,[sol_cut])
                append!(solution_cut_u,[sol_cut_u])
            end
        end    
        
        for l in 1:length(ax)-1
            if length(ax[l+1][:lines])>0 #clear plot in current axis in case there is any
                ax[l+1][:lines] = []
            end
        end
        
        if cut_type=="solutions" 
            for l in 1:nvars      
                ax[l+1].plot(Z, remove_singleton!(solution_cut_s[l,:,:]'),lw=3)
                ax[l+1].plot(Z, remove_singleton!(solution_cut_u[l,:,:]'),ls="--",lw=3)
            end    
        elseif cut_type=="jacobian_eigenvalues"
            for l in 1:nsolsmax
                ax[l+1].axhline(0,ls=":",lw=4) #reference value for unstability test 
                for sols in remove_singleton!.([solution_cut_s[:,l,:]',solution_cut_u[:,l,:]'])
                    _plot_2D_solutions_jacobian_cut(res,sols,parameter_val,Z,ax[l+1])         
                end
            end
        elseif cut_type=="transform"  
            for l in 1:length(string_f)
                ax[l+1].plot(Z,remove_singleton!(solution_cut_s[l]'),ls="-",lw=3)
                ax[l+1].plot(Z, remove_singleton!(solution_cut_u[l]'),ls="--",lw=3)
            end
        end    
        
        legend_elements = [plt.Line2D([1], [1], linestyle="-", color="k", label=lab[1],
                        markerfacecolor="k", markersize=5),
                        plt.Line2D([1], [1], linestyle="--", color="k", label=lab[2],
                        markerfacecolor="k", markersize=5)]    
        ax[end].legend(handles=legend_elements,loc="best") 
    end    

    f[:canvas][:mpl_connect]("button_press_event", onclick)    
    
    #When the interactive plotting window is closed, subsequent plots are not interactive if not requested so
    function on_close(event)
        pygui(false) 
    end
    f[:canvas][:mpl_connect]("close_event",on_close)

    ax = resize_axes!(f,ax,nrows,ncols)
    if nrows>1
        f.colorbar(im, ax=ax) #common colorbar if a list of axes is passed to colorbar instead of a single one
    else
        f.colorbar(im, ax=ax[end])
    end
end

############################################################
#DEPRECATED
#"Create latex labels for subsequent plots"
#function get_latex_label(x::String)
#    if REPL.symbol_latex(x) == "" #the character is not a symbol
#        latex_label = string("\$",x,"\$")
#    else
#        latex_label = string("\$",REPL.symbol_latex(x),"\$")
#    end
#    latex_label
#end


#= function mod_phase_plots(phase_diagram::phase_diagram)
    """Plot abs,angle for all physical_solutions, including unstable ones"""
    nvar  = size(phase_diagram.physical_solutions)[1] #number of variables
    nsols = size(phase_diagram.physical_solutions)[2]
    sweep1 = reduce(hcat,collect(values(phase_diagram.p_sweep)))[:,1] #values for the swept parameter
    sweep2 = reduce(hcat,collect(values(phase_diagram.p_sweep)))[:,2] #values for the swept parameter
    f,axx = subplots(nsols,nvar,figsize=(4*nvar,4*nsols))
    
    Us   = [view(phase_diagram.physical_solutions, m, :,:,:)   for m in 1:2:nvar] #view prevents the copy of the array
    Vs   = [view(phase_diagram.physical_solutions, m, :,:,:)   for m in 2:2:nvar]
    Xs   = [sqrt.(us.^2 + vs.^2) for (us,vs) in zip(Us,Vs)]
    Args = [atan.(vs./us) for (us,vs) in zip(Us,Vs)]
    
    if nsols==1
        axx = add_dim(axx)
    end
    print(length(axx))

    extent = [sweep1[1],sweep1[end],sweep2[1],sweep2[end]]
    for k in 1:Int(floor(nvar/2))
        for l in 1:nsols
            axx[l,k].imshow(Xs[k][l,:,:],extent=extent,aspect="auto",vmin=0,vmax=2)
            axx[l,k+Int(floor(nvar/2))].imshow(Args[k][l,:,:],extent=extent,aspect="auto",cmap="hsv",vmin=0,vmax=1)
            axx[l,k].set_xlabel(string("\$",collect(keys(phase_diagram.p_sweep))[1],"\$"),fontsize=24)
            axx[l,k].set_ylabel(string("\$",collect(keys(phase_diagram.p_sweep))[2],"\$"),fontsize=24)
            axx[l,k+Int(floor(nvar/2))].set_xlabel(string("\$",collect(keys(phase_diagram.p_sweep))[1],"\$"),fontsize=24)
            axx[l,k+Int(floor(nvar/2))].set_ylabel(string("\$",collect(keys(phase_diagram.p_sweep))[2],"\$"),fontsize=24)
        end
        axx[1,k].set_title(string("\$X_",string(k),"\$"))
        axx[1,k+Int(floor(nvar/2))].set_title(string(L"$\Phi_",string(k),"\$"))
    end
    f.tight_layout()  
end =#
