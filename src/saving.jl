"""
    $(SIGNATURES)

Saves `object` into `.jld2` file `filename` (the suffix is added automatically if not entered).
The resulting file contains a dictionary with a single entry.

"""
function save(filename, object)
    return JLD2.save(_jld2_name(filename), Dict("object" => object))
end

function save(filename, x::Result)
    x_nofunc = deepcopy(x)

    # compiled functions cause problems in saving: ignore J now, compile when loading
    x_nofunc.jacobian = 0
    return JLD2.save(_jld2_name(filename), Dict("object" => x_nofunc))
end

_jld2_name(filename) = filename[(end - 4):end] == ".jld2" ? filename : filename * ".jld2"

"""
    $(SIGNATURES)

Loads an object from `filename`. For objects containing symbolic expressions such as `HarmonicEquation`, the symbolic variables are
reinstated in the `HarmonicBalance` namespace.

"""
function load(filename)
    loaded = JLD2.load(filename)
    if haskey(loaded, "object") #otherwise save data is from a plot
        loaded = loaded["object"]

        # we need the symbols in our namespace to parse strings with `transform_solutions`
        _parse_symbol_names(loaded)

        # automatic saving fails for some objects: reconstruct these manually
        _parse_loaded(loaded)
    else
        return loaded
    end
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
    all_symbols = cat(
        x.parameters,
        get_variables(x.variables),
        get_independent_variables(x.eom),
        get_independent_variables(x.eom.natural_equation);
        dims=1,
    )
    return declare_variable.(string.(all_symbols))
end

_parse_symbol_names(x::Result) = _parse_symbol_names(x.problem)
_parse_symbol_names(x) = nothing

# Exporting to csv

"""
    $(SIGNATURES)

Saves into `filename` a specified solution `branch` of the Result `res`.
"""
function export_csv(filename::String, res::Result, branch::Int)
    branch_data = getindex.(res.solutions, branch)
    return writedlm(filename, branch_data, ',')
end

export_csv(filename::String, res::Result; branch::Int64) = export_csv(filename, res, branch)
