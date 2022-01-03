using PyPlot
using PyCall
using Latexify
export plot_2D_phase_diagram_interactive

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
function _get_interactive_plot_axes(x,y,gx,gy,var_names,cut_dim,cut_type,nvars,nsolsmax_physical,sol_type,not_sol_type; string_f)
    if  cut_type=="solutions"    
        N_panels = nvars + 1
        lab = [sol_type,not_sol_type]
    elseif cut_type=="jacobian_eigenvalues"
        N_panels = nsolsmax_physical + 1
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
    
    return sc,ax,f,annot,lab
end


"Plotter of Jacobian eigenvalues along a cut of a 2D solutions. TODO: condense this with  plot_1D_jacobian_eigenvalues"
function _plot_2D_solutions_jacobian_cut(res::Result,filtered_sol,parameter_cut,Z::Vector{Float64},ax::PyCall.PyObject)
    
     # pre-evaluate the Jacobian with fixed parameters
    fixed_pair     = parameter_cut[2]
    fixed_params = Dict(k=>parse(ComplexF64,string(v))  for (k,v) in pairs(merge(res.fixed_parameters,Dict(fixed_pair))))
    fixed_params = Num_to_Variable(fixed_params)
    
    Jac = HomotopyContinuation.evaluate(HomotopyContinuation.jacobian(res.problem.system),collect(keys(fixed_params))=>collect(values(fixed_params))) #TODO make this work with HarmonicBalance.Jacobian
    
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

    nsolsmax_physical = sum(any.(classify_branch(res, "physical"))) #number of physical solutions
    sc,ax,f,annot,lab = _get_interactive_plot_axes(x,y,gx,gy,var_names,cut_dim,cut_type,nvars,nsolsmax_physical,sol_type,not_sol_type,string_f=string_f)
    
    length(vec(ax)) <= nrows*ncols || error("insufficient # of panels requested, please increase nrows or ncols") #sanity check before any plot is made
   
    im,Nmax = plot_2D_phase_diagram(res; stable=stable,observable=observable,ax=ax[1])

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
        _prepare_colorbar(f,ax,im,Nmax)
    else
        _prepare_colorbar(f,ax[end],im,Nmax)
    end
end