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
function ExprUtils.substitute_all(eq::Union{Num,Equation}, rules::Dict{HarmonicVariable})
    return Symbolics.substitute(
        eq, Dict(zip(getfield.(keys(rules), :symbol), values(rules)))
    )
end

function ExprUtils.substitute_all(var::HarmonicVariable, rules)
    sym, freq = var.symbol, var.ω
    return HarmonicVariable(
        substitute_all(sym, rules),
        var.name,
        var.type,
        substitute_all(freq, rules),
        var.natural_variable,
    )
end

function ExprUtils.substitute_all(vars::Vector{HarmonicVariable}, rules)
    return [substitute_all(var, rules) for var in vars]
end

"Returns the symbols of a `HarmonicVariable`."
get_variables_nums(vars::Vector{Num}) =
    unique(flatten([Num.(get_variables(x)) for x in vars]))

Symbolics.get_variables(var::HarmonicVariable)::Num = Num(first(get_variables(var.symbol)))

Base.isequal(v1::HarmonicVariable, v2::HarmonicVariable)::Bool =
    isequal(v1.symbol, v2.symbol)

"The derivative of f w.r.t. x of degree deg"
function d(f::Num, x::Num, deg=1)::Num
    return isequal(deg, 0) ? f : (Differential(x)^deg)(f)
end
d(funcs::Vector{Num}, x::Num, deg=1) = Num[d(f, x, deg) for f in funcs]

"Declare a variable in the the currently active namespace"
function declare_variable(name::String)
    var_sym = Symbol(name)
    @eval($(var_sym) = first(Symbolics.@variables $var_sym))
    return eval(var_sym)
end

"Declare a variable that is a function of another variable in the the current namespace"
function declare_variable(name::String, independent_variable::Num)
    # independent_variable = declare_variable(independent_variable) convert string into Num
    var_sym = Symbol(name)
    new_var = Symbolics.@variables $var_sym(independent_variable)
    @eval($(var_sym) = first($new_var)) # store the variable under "name" in this namespace
    return eval(var_sym)
end

"Return the name of a variable (excluding independent variables)"
function var_name(x::Num)
    var = Symbolics._toexpr(x)
    return var isa Expr ? String(var.args[1]) : String(var)
end
#  var_name(x::Term) = String(Symbolics._toexpr(x).args[1])
var_name(x::SymbolicUtils.Sym) = String(x.name)
