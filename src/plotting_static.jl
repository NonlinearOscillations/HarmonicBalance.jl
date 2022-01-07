using PyPlot
using PyCall
using Latexify
using JLD2
export plot_1D_solutions, plot_2D_phase_diagram, transform_solutions
export _set_plotting_settings, _prepare_colorbar, _prettify_label


"Set global plotting settings"
function _set_plotting_settings()
    plt.style.use("default") #reset settings
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams") 
    rcParams["text.usetex"]       = true
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
    rcParams["figure.dpi"] = 220
end


# the total number of solutions
_total_elements(array) = prod(size(array))

_add_dim!(x::Array) = reshape(x, (1,size(x)...)) #quickfix for the case where axx is not a matrix. It would require an extra singleton dimension

_squeeze!(arr) = dropdims(arr, dims = tuple(findall(size(arr) .== 1)...)) #remove singleton dimensions of an array

"""
$(TYPEDSIGNATURES)

Goes over a solution and an equally-sized array (a "mask") of booleans. 
true  -> solution unchanged
false -> changed to NaN (omitted from plotting)
"""
function filter_solutions(solution::Vector,  booleans)
    _total_elements(solution) == _total_elements(booleans) || error("attempt to filter a solution using a wrongly-sized boolean array")
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

"""
$(TYPEDSIGNATURES)

For a transformed, physical solution, evaluates multi-solution condition (e.g. Base.maximum or Base.argmax).
If real_function=true (default) the function is assume to be real and small imaginary parts are discarded.
"""
function map_multi_solutions(solution,mapping=nothing; real_function=true)
    length(vcat(solution...)[1])>1  && error("input solutions have wrong dimensions. Please transform first.")#sanity check for solutions still containing the solution tuple u1,v1,... (not transformed) 
    isempty(methods(mapping)) && error("mapping is not a Julia function") #sanity check for mapping being a function
    if real_function
        @warn "imaginary parts discarded in multi-solution map"
        return mapping.(real.(solution))
    else
        return mapping.(solution)
    end
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

"inserts underscores in variable names for prettier transformed laTeX strings"
function _prettify_label(res::Result,label::String)
    replace_rules = [string("u",k)=>string("u_",k) for k in 1:length(res.problem.variables)÷2]
    append!(replace_rules,[string("v",k)=>string("v_",k) for k in 1:length(res.problem.variables)÷2])
    return reduce(replace, replace_rules, init=label)
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

    f,ax = subplots(1,1,figsize=(7,4))

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

    
    if !isnothing(filename)
        xdata,ydata = [line.get_xdata() for line in lines], [line.get_ydata() for line in lines]
        markers = [line.get_marker() for line in lines]
        marker_dict = Dict(marker=>sol_type,"X"=>not_sol_type)
        JLD2.save(_jld2_name(filename), Dict(string(x) => xdata,string(y)=>ydata,"marker_dict"=>marker_dict,"markers"=>markers))
    end

    ax.set_xlabel(Latexify.latexify(x),fontsize=24) 
    ax.set_ylabel(Latexify.latexify(_prettify_label(res,y)),fontsize=24) 

    #legend preparation
    ignored_idx = [all(isnan.(line.get_ydata())) for line in lines] #make up a legend with only non ignored entries in the plotter
    Nb   = sum(.~ignored_idx) #number of branches
    leg1 = ax.legend(string.(collect(1:Nb)),ncol=(Nb<10) + Nb÷10*(Nb>10),
            bbox_to_anchor=(Nb÷10 + 1, 0.95))
    ax.add_artist(leg1)
    
    h_leg2 = [plt.Line2D([0], [0], marker=marker, color="w", label=sol_type,    markerfacecolor="k", markersize=10),
            plt.Line2D([0], [0], marker="X"   , color="w", label=not_sol_type,markerfacecolor="k", markersize=10)] 

    ax.legend(handles=h_leg2,bbox_to_anchor=(-0.25, 0.95)) 
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    f.tight_layout()
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
- `ax`: axis object from `PyCall.PyObject` setting the coordinate system where data will be plotted. If not given, it is created automatically.
- `filename`: if different from `nothing`, plotted data and parameter values are exported to `./filename.jld2`.

"""
function plot_1D_jacobian_eigenvalues(res::Result; x::String, physical=true, stable=false,marker_re="o",marker_im="X",ax=nothing, filename=nothing)
    _set_plotting_settings()

    xplot = transform_solutions(res, x) #independent variable to plot
    
    #filtering of the solutions according to keyword arguments
    !physical && stable && error("Stability is not defined for unphysical solutions!")
    
    numbers = [Base.not_int.(isnan.(prod.(s))) for s in res.solutions] # true for solution which are not NaN
    to_evaluate = physical ? (stable ? res.classes["stable"] : res.classes["physical"]) : numbers
    to_evaluate = [to_evaluate[i] .* numbers[i] for (i,_) in enumerate(numbers)] 


    nsolsmax  = sum(any.(classify_branch(res, "physical"))) #maximum number of physical solutions
    input_ax = ax
    if isnothing(input_ax) #create figure if axes are not provided, otherwise accept input
        f,ax = subplots(1,nsolsmax,figsize=(4*nsolsmax,4))
    end

    if nsolsmax==1
        ax = _add_dim!([ax])
    end

    lines_re,lines_im =[],[]
    for branch in 1:nsolsmax
        ax[branch].axhline(0,ls="dashdot",lw=2,alpha=0.8) #reference value for unstability test 
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
   ax[1].legend(handles=legend_elements,loc="best",fontsize=15) 
end

"""
    plot_2D_solutions(res::Result,ax=nothing; filename=nothing)

Make a 2D plot of each of solutions vs swept parameters for a `Result` object, obtained from a `get_steady_states` applied to a 2D parameter grid.`.

Keyword arguments

- `ax`: axis object from `PyCall.PyObject` setting the coordinate system where data will be plotted. If not given, it is created automatically.
- `filename`: if different from `nothing`, plotted data and parameter values are exported to `./filename.jld2`.
"""
function plot_2D_solutions(res::Result; ax=nothing, filename=nothing, z=nothing,plot_only="nothing")
    _set_plotting_settings()
    nvar  = length(res.solutions[1,1][1]) #number of variables
    nsols = length(res.solutions[1,1]) #maximum number of solutions
    x,y = collect(values(res.swept_parameters))
    extent=[x[1],x[end],y[1],y[end]]

    if isnothing(z)
        Z = res.solutions
    else
        Z =   transform_solutions(res, z) # first transform, then filter
    end
    #transform the solution structs into tensors that are easier to handle by imshow
    physical_solutions = real.(filter_solutions.(Z,res.classes["physical"]))
    if plot_only=="stable"
        physical_solutions = real.(filter_solutions.(Z,res.classes["stable"]))
    end

    if isnothing(z)
        physical_sols = reshape(reduce(hcat,[reduce(hcat,sol) for sol in physical_solutions]),nvar,nsols,length(x),length(y))
    else
        physical_sols = reshape(reduce(hcat,[reduce(hcat,sol) for sol in physical_solutions]),nsols,length(x),length(y))
    end
    var_names = _get_var_name_labels(res) #variable names for plot labels

    input_ax = ax
    if isnothing(input_ax) #create figure if axes are not provided, otherwise accept input
        if isnothing(z)
            f,ax = subplots(nsols,nvar,figsize=(4*nvar,4*nsols))
        else
            f,ax = subplots(1,nsols,figsize=(4*nsols,4))
        end    
    end

    px,py = string.(collect(keys(res.swept_parameters)))  #swept parameter strings
    if isnothing(z)
        save_dict = Dict([string("panel (",m,",",l,")")=> Dict() for m in 1:nvar for l in 1:nsols])
        for m in 1:nvar
            for l in 1:nsols 
                a = ax[l,m].imshow(physical_sols[m,l,:,end:-1:1]',extent=extent,aspect="auto")
                colorbar(a,ax=ax[l,m])
                ax[l,m].set_xlabel(Latexify.latexify(px),fontsize=24); 
                ax[l,m].set_ylabel(Latexify.latexify(py),fontsize=24); 
                if !isnothing(filename)
                    save_dict[string("panel (",m,",",l,")")]= Dict("variable"=>var_names[m],"solution #"=>l,"data"=>a.get_array(),
                    string("(",px,"_min ",px,"_max ",py,"_min ",py,"_max)")=>extent)
                end
            end
            ax[1,m].set_title(Latexify.latexify(_prettify_label(res,var_names[m])),fontsize=20)
        end
    else
        save_dict = Dict([string("panel (",l,")")=> Dict() for l in 1:nsols])
        for l in 1:nsols 
            a = ax[l].imshow(physical_sols[l,:,end:-1:1]',extent=extent,aspect="auto")
            colorbar(a,ax=ax[l])
            ax[l].set_xlabel(Latexify.latexify(px),fontsize=24); 
            ax[l].set_ylabel(Latexify.latexify(py),fontsize=24); 
            if !isnothing(filename)
                save_dict[string("panel (",l,")")]= Dict("variable"=>z,"solution #"=>l,"data"=>a.get_array(),
                string("(",px,"_min ",px,"_max ",py,"_min ",py,"_max)")=>extent)
            end
            ax[l].set_title(string("solution ",l),fontsize=18)
        end
        f.suptitle(Latexify.latexify(_prettify_label(res,z)),fontsize=20,y=1.01)
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
    Nmax = maximum(obs_2D)
    im = ax.imshow(obs_2D[:,end:-1:1]', extent=extent, aspect="auto",vmin=minimum(obs_2D),vmax=Nmax)
    if isnothing(input_ax) 
        _prepare_colorbar(f,ax,im,Nmax)
    end

    px,py =string.([x,y])#swept parameter strings

    ax.set_xlabel(Latexify.latexify(px),fontsize=24)
    ax.set_ylabel(Latexify.latexify(py),fontsize=24)
    ax.set_title(pd_title)

    if !isnothing(filename)
        JLD2.save(_jld2_name(filename), Dict("observable"=>observable,"data"=>im.get_array(),
                                                 string("(",px,"_min ",px,"_max ",py,"_min ",py,"_max)")=>extent))
    end
    im,Nmax
end

function _prepare_colorbar(f,ax,im,Nmax;Nmax_discrete=9)  #creates and inserts in a phase diagram figure f a discrete or continuous colorbar
    if Nmax <= Nmax_discrete
        f.colorbar(im, ax=ax,ticks=collect(1:Nmax), boundaries=collect(1:Nmax+1).-0.5)
    else
        f.colorbar(im, ax=ax)
    end
end


