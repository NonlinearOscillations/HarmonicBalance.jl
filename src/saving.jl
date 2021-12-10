using DataStructures
using JLD2


function save(filename, x)
    JLD2.save(filename, Dict("object" => x))
end


function save(filename, x::Result)
    x_nofunc = deepcopy(x)
    # compiled functions cause problems in saving: ignore J now, compile when loading
    x_nofunc.jacobian = false
    JLD2.save(filename, Dict("object" => x_nofunc))
end

function load(filename)
    loaded = JLD2.load(filename)["object"]
    
    # we need the symbols in our namespace to parse strings with `transform_solutions`
    _parse_symbol_names(loaded)

    # automatic saving fails for some objects: reconstruct these manually
    _parse_loaded(loaded)
end


function _parse_loaded(x::Problem)
    # reconstruct the HomotopyContinuation System
    system = HC_wrapper.System(x.eom)
    x.system = system
    return x
end


function _parse_loaded(x::Result)
    # reconstruct System and the compiled Jacobian
    x.problem.system = HC_wrapper.System(x.problem.eom)
    x.jacobian = _compile_Jacobian(x.problem, x.swept_parameters, x.fixed_parameters)
    return x
end

_parse_loaded(x) = x


"Retrieve names for all symbols and declare them in this namespace."
function _parse_symbol_names(x::Problem)
    all_symbols= cat(x.parameters, get_variables(x.variables), 
    get_independent_variables(x.eom), get_independent_variables(x.eom.natural_equation), dims=1)
    return declare_variable.(string.(all_symbols))
end

_parse_symbol_names(x::Result) = _parse_symbol_names(x.problem)
_parse_symbol_names(x) = nothing
#=


"Retrieve names for all symbols"
function dump_symbol_names(r::Result)
    all_symbols= cat(r.problem.parameters, get_variables(r.problem.variables), dims=1)
    return string.(all_symbols)
end


"Returns an unique name for a given solution object, exploiting DrWatson library"
function get_filename(fixed_parameters,swept_parameters)
    name_dict = copy(fixed_parameters) #this could be potentially simplified by the macro @strdict, which creates a dictionary from existing Variable
    for (parameter,value) in zip(keys(swept_parameters),values(swept_parameters))
        name_dict[declare_variable(string(parameter,"_min"))] = value[1]
        name_dict[declare_variable(string(parameter,"_max"))] = value[end]
    end
    return savename(name_dict,"jld2")
end



"""
    save_result(project_name,swept_parameters::ParameterRange, fixed_parameters::ParameterList,res::Result; filename=nothing)

Saves a `Result` struct for a steadystate problem, together with the parameters of the simulation into a .jld2 file stored in the folder `./project_name`. 
This format comprises a subset of `HDF5`, see `JLD2.jl` documentation.    
If `filename` is not specified, it is determined automatically from the simulation parameters. 
If the filename already exists, an incremental backup-number is introduced (e.g. running as #1,#2,#3...). 
See documentation of `safesave` function in `DrWatson.jl` for further details.`

This function is called in `get_steady_states`.
"""
function save_result(project_name::String,swept_parameters::ParameterRange, fixed_parameters::ParameterList,res::Result; filename=nothing)
    if isnothing(filename) 
        filename_str = get_filename(fixed_parameters,swept_parameters)
    else
        filename_str = string(filename,".jld2")
    end
    swept_params_str = OrderedDict(zip(string.(keys(swept_parameters)),values(swept_parameters)))
    fixed_params_str = OrderedDict(zip(string.(keys(fixed_parameters)),values(fixed_parameters)))

    solutions, classes = res.solutions, res.classes
    names = dump_symbol_names(res)
    savefile = @strdict solutions classes swept_params_str fixed_params_str names #create saving struct with only standard types (nothing symbolic)

    safesave(datadir(project_name, filename_str), savefile) #if overwritting is allowed, then wsave is enough
end


"""
    load_result(project_name::String,index::Int64, prob::Problem)

Lists all saved data corresponding to a known `Problem`, stored  at `./project_name`. 
Instantiates the `Result` class with the solution data, parameter values and stability analyses arrays, stored at the position `index` in the path list.  

Example: loading data from the last executed simulation, stored in `./test_folder` 
```julia
project_name = "test_folder" #used to create data folders later on
idx_sim = 1 #simulation index
load_sol = HarmonicBalance.load_result(project_name,idx_sim,EOM_HC)
```
"""
function load_result(project_name::String,index::Int64, prob::Problem; fixed_param_input=nothing)
    simulation_path = readdir(datadir(project_name))[index] 
    simulation_data = wload(datadir(project_name,simulation_path))
    declare_variable.(simulation_data["names"])
    display(simulation_data)    
    swept_params    = ParameterRange(zip(declare_variable.(keys(simulation_data["swept_params_str"])),values(simulation_data["swept_params_str"]))) #sympy Variable in swept_params_dictionary need to be stored as strings. This is the reverse transformation to Variable
    if haskey(simulation_data, "fixed_parameters_str")
         fixed_params    = ParameterList(zip(declare_variable.(keys(simulation_data["fixed_params_str"])),values(simulation_data["fixed_params_str"])))
    else
        fixed_params = fixed_param_input
    end
    Result(simulation_data["solutions"], swept_params, fixed_params, prob, simulation_data["classes"])
end
=#