# pretty-printing
Base.display(var::HarmonicVariable) = display(var.name)
Base.display(var::Vector{HarmonicVariable}) = display.(getfield.(var, Symbol("name")))

function _coordinate_transform(new_var, ω, t, type)::Num
    coords = Dict([
        "u" => new_var * cos(ω * t), "v" => new_var * sin(ω * t), "a" => new_var
    ])
    return coords[type]
end

function _create_harmonic_variable(
    nat_var::Num, ω::Num, t::Num, type::String; new_symbol::String
)::Tuple{Num,HarmonicVariable}
    new_var = declare_variable(new_symbol, t) # this holds the internal symbol
    name = type * "_{" * var_name(nat_var) * "," * Base.replace(string(ω), "*" => "") * "}"
    rule = _coordinate_transform(new_var, ω, t, type) # contribution of this harmonic variable to the natural variable
    hvar = HarmonicVariable(new_var, name, type, ω, nat_var)
    return rule, hvar
end

###
# Functions for variable substutions and manipulation of HarmonicVariable
###

# when HV is used for substitute, substitute its symbol
function substitute_all(eq::Union{Num,Equation}, rules::Dict{HarmonicVariable})
    return substitute(eq, Dict(zip(getfield.(keys(rules), :symbol), values(rules))))
end

function substitute_all(var::HarmonicVariable, rules)
    sym, freq = var.symbol, var.ω
    return HarmonicVariable(
        substitute_all(sym, rules),
        var.name,
        var.type,
        substitute_all(freq, rules),
        var.natural_variable,
    )
end

function substitute_all(vars::Vector{HarmonicVariable}, rules)
    return [substitute_all(var, rules) for var in vars]
end

"Returns the symbols of a `HarmonicVariable`."
get_variables_nums(vars::Vector{Num}) =
    unique(flatten([Num.(get_variables(x)) for x in vars]))

Symbolics.get_variables(var::HarmonicVariable)::Num = Num(first(get_variables(var.symbol)))

Base.isequal(v1::HarmonicVariable, v2::HarmonicVariable)::Bool =
    isequal(v1.symbol, v2.symbol)
