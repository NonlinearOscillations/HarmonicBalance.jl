import Symbolics.get_variables; export get_variables

# pretty-printing
display(vars::HarmonicVariable) = display(vars.names)
display(vars::Vector{HarmonicVariable}) = display.(getfield.(vars, Symbol("names")))


# Each keyword gives two names for the new variables
coordinate_names = Dict([["Cartesian", ["u", "v"]], ["polar", ["a", "ϕ"]]])


function _coordinate_transform(new_vars, ω, t)
    Dict([["Cartesian",  new_vars[1] * cos(ω*t) + new_vars[2] * sin(ω*t)],
["polar", new_vars[1] * cos(ω*t + new_vars[2])]])
end


"Generates two variables describing the harmonic ω of a natural variable."
function _rotate_variable(nat_var::Num, ω::Num, time::Num, transform::String; new_symbols::Vector{String})

    types = coordinate_names[transform] # is the transformation "Cartesian" or "polar"
    new_vars = [declare_variable(s, time) for s in new_symbols] # this holds the internal symbols

    var_names_strings = [VDP_name * "_{" * var_name(nat_var) * "," * Base.replace(string(ω),"*"=>"") * "}" for VDP_name in types]
    var_names = Dict(zip([Num(s) for s in new_symbols], var_names_strings))
    repl_rule = _coordinate_transform(new_vars, ω, time)[transform] # a dictionary mapping the coordinate type to an expression
    VDP = HarmonicVariable(new_vars, var_names, types, ω, nat_var)

    repl_rule, VDP
end


###
# Functions for variable substutions and manipulation of HarmonicVariable
###

function substitute_all(vars::HarmonicVariable, rules)
    sym, freq = vars.symbols, vars.ω
    HarmonicVariable(substitute_all(sym, rules), vars.names, vars.types, substitute_all(freq, rules), vars.natural_variable)
end


substitute_all(vars::Vector{HarmonicVariable}, rules) = [substitute_all(var, rules) for var in vars]


"Returns the symbols of a `HarmonicVariable`."
get_variables(vars::Vector{Num}) = unique(flatten([Num.(get_variables(x)) for x in vars]))

get_variables(vars::HarmonicVariable) = unique(get_variables(vars.symbols))





