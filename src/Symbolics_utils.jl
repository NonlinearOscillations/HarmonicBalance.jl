using Symbolics.SymbolicUtils: add_with_div, frac_similarterm

export rearrange
export drop_powers
export get_averaged_equations
export d
export substitute_all
export get_all_terms
export var_name


"The derivative of f w.r.t. x of degree deg"
function d(f::Num, x::Num, deg=1)
    isequal(deg,0) ? f : (Differential(x)^deg)(f)
end
d(funcs::Vector{Num}, x::Num, deg=1) = [d(f, x, deg) for f in funcs]

"Declare a variable in the the currently active namespace"
 function declare_variable(name::String)
    var_sym = Symbol(name)
    @eval($(var_sym) = first(@variables $var_sym))
    return eval(var_sym)
 end


 "Declare a variable that is a function of another variable in the the current namespace"
 function declare_variable(name::String, independent_variable::Num)
    # independent_variable = declare_variable(independent_variable) convert string into Num
    var_sym = Symbol(name)
    new_var = @variables $var_sym(independent_variable)
    @eval($(var_sym) = first($new_var)) # store the variable under "name" in this namespace
    return eval(var_sym)
 end


 "Return the name of a variable (excluding independent variables)"
function var_name(x::Num)
    var = Symbolics._toexpr(x)
    return var isa Expr ? String(var.args[1]) : String(var)
end
#  var_name(x::Term) = String(Symbolics._toexpr(x).args[1])
 var_name(x::Sym) = String(x.name)


 """
 $(TYPEDSIGNATURES)

 Perform substitutions in `rules` on `x`.
 `include_derivatives=true` also includes all derivatives of the variables of the keys of `rules`.
 """
 function substitute_all(x::T, rules::Dict; include_derivatives=true)::T where {T<:Union{Equation, Num}}
    if include_derivatives
        rules = merge(rules, Dict([Differential(var) => Differential(rules[var]) for var in keys(rules)]))
    end
    return substitute(x, rules)
 end

 "Variable substitution - dictionary"
function substitute_all(dict::Dict, rules::Dict)::Dict
    new_keys = substitute_all.(keys(dict), rules)
    new_values = substitute_all.(values(dict), rules)
    return Dict(zip(new_keys, new_values))
end

substitute_all(v::Union{Array{Num}, Array{Equation}}, rules::Union{Dict, Pair, Vector}) = [substitute_all(x, rules) for x in v]
substitute_all(x::Union{Num, Equation}, rules::Union{Pair, Vector, Dict}) = substitute_all(x, Dict(rules))
substitute_all(x, rules::OrderedDict) = substitute_all(x, Dict(rules))
substitute_all(x::Complex{Num}, rules) = substitute_all(Num(x.re.val.arguments[1]), rules)
substitute_all(x::Symbolics.Arr, args...) = [substitute_all[r] for i in x]
substitute_all(x, rules) = substitute_all(Num(x), rules)


"""
$(SIGNATURES)
Remove parts of `expr` where the combined power of `vars` is => `deg`.

# Example
```julia-repl
julia> @variables x,y;
julia>drop_powers((x+y)^2, x, 2)
y^2 + 2*x*y
julia>drop_powers((x+y)^2, [x,y], 2)
0
julia>drop_powers((x+y)^2 + (x+y)^3, [x,y], 3)
x^2 + y^2 + 2*x*y
```
"""
function drop_powers(expr::Num, vars::Vector{Num}, deg::Int)
    @variables ϵ
    subs_expr = deepcopy(expr)
    rules = Dict([var => ϵ*var for var in unique(vars)])
    subs_expr = Symbolics.expand(substitute_all(subs_expr, rules))
    max_deg = max_power(subs_expr, ϵ)
    removal = Dict([ϵ^d =>Num(0) for d in deg:max_deg])
    res = substitute_all(substitute_all(subs_expr, removal), Dict(ϵ => Num(1)))
    Symbolics.expand(res)
    #res isa Complex ? Num(res.re.val.arguments[1]) : res
end


drop_powers(expr::Vector{Num}, var::Num, deg::Int) = [drop_powers(x, var, deg) for x in expr]

# calls the above for various types of the first argument
drop_powers(eq::Equation, var, deg) = drop_powers(eq.lhs, var, deg) .~ drop_powers(eq.lhs, var, deg)
drop_powers(eqs::Vector{Equation}, var, deg) = [Equation(drop_powers(eq.lhs, var, deg), drop_powers(eq.rhs, var, deg)) for eq in eqs]
drop_powers(expr, var::Num, deg::Int) = drop_powers(expr, [var], deg)
drop_powers(x, vars, deg) = drop_powers(Num(x), vars, deg)


flatten(a) = collect(Iterators.flatten(a))


###
# STUFF BELOW IS MAINLY FOR FOURIER-TRANSFORMING
###


get_independent(x::Num, t::Num) = get_independent(x.val, t)
get_independent(x::Complex{Num}, t::Num) = get_independent(x.re, t) + im*get_independent(x.im, t)
get_independent(v::Vector{Num}, t::Num) = [get_independent(el, t) for el in v]
get_independent(x, t::Num) = x

function get_independent(x::BasicSymbolic, t::Num)
    if isadd(x)
        return sum([get_independent(arg, t) for arg in arguments(x)])
    elseif ismul(x)
        return prod([get_independent(arg, t) for arg in arguments(x)])
    elseif ispow(x)
        return !is_function(x.base, t) && !is_function(x.exp, t) ? x : 0
    elseif isdiv(x)
        return !is_function(x.den, t) ? get_independent(x.num, t) / x.den : 0
    elseif isterm(x) || issym(x)
        return !is_function(x, t) ? x : 0
    else
        return x
    end
end

"Return all the terms contained in `x`"
get_all_terms(x::Num) = unique(_get_all_terms(Symbolics.expand(x).val))
get_all_terms(x::Equation) = unique(cat(get_all_terms(Num(x.lhs)), get_all_terms(Num(x.rhs)), dims=1))

_get_all_terms_mul(x) = Num.(SymbolicUtils.arguments(x))
_get_all_terms_div(x) = Num.([_get_all_terms(x.num)..., _get_all_terms(x.den)...])
_get_all_terms(x) = Num(x)

function _get_all_terms_add(x)::Vector{Num}
    list = []
    for term in keys(x.dict)
        list = cat(list, _get_all_terms(term), dims=1)
    end
    return list
end

function _get_all_terms(x::BasicSymbolic)
    if isadd(x)
        return _get_all_terms_add(x)
    elseif ismul(x)
        return _get_all_terms_mul(x)
    elseif isdiv(x)
        return _get_all_terms_div(x)
    else
        return Num(x)
    end
end

function is_harmonic(x::Num, t::Num)::Bool
    all_terms = get_all_terms(x)
    t_terms = setdiff(all_terms, get_independent(all_terms, t))
    isempty(t_terms) && return true
    trigs = is_trig.(t_terms)

    if !prod(trigs)
        return false
    else
        powers = [max_power(first(term.val.arguments), t) for term in t_terms[trigs]]
        return all(isone, powers)
    end
end

is_harmonic(x::Equation, t::Num) = is_harmonic(x.lhs, t) && is_harmonic(x.rhs, t)
is_harmonic(x, t) = is_harmonic(Num(x), Num(t))

"Convert all sin/cos terms in `x` into exponentials."
function trig_to_exp(x::Num)
    all_terms = get_all_terms(x)
    trigs = filter(z -> is_trig(z), all_terms)

    rules = []
    for trig in trigs
        is_pow = ispow(trig.val) # trig is either a trig or a power of a trig
        power = is_pow ? trig.val.exp : 1
        arg = is_pow ? arguments(trig.val.base)[1] : arguments(trig.val)[1]
        type = is_pow ? operation(trig.val.base) : operation(trig.val)

        if type == cos
            term = Complex{Num}((exp(im*arg) + exp(-im*arg))^power * (1//2)^power,0)
        elseif type == sin
            term = (1*im^power)* Complex{Num}( ((exp(-im*arg) - exp(im*arg)))^power*(1//2)^power, 0)
        end
        # avoid Complex{Num} where possible as this causes bugs
        # instead, the Nums store SymbolicUtils Complex types
        term = Num(Symbolics.expand(term.re.val + im*term.im.val))
        append!(rules, [trig => term])
    end

    result = Symbolics.substitute(x, Dict(rules))
    #result = result isa Complex ? Num(first(result.re.val.arguments)) : result
    result = Num(result)
    return result
end


"Return true if `f` is a function of `var`."
is_function(f, var) = any(isequal.(get_variables(f), var))


"Return true if `f` is a sin or cos."
function is_trig(f::Num)
    f = ispow(f.val) ? f.val.base : f.val
    isterm(f) && SymbolicUtils.operation(f) ∈ [cos,sin] && return true
    return false
end


"A vector of Sym(0) of length n"
Num_zeros(n::Int64) = [Num(0) for k in 1:n]
Num_zeros(vec::Vector{Any}) = Num_zeros(length(vec))


"""
$(TYPEDSIGNATURES)
Returns the coefficient of cos(ωt) in `x`.
"""
function fourier_cos_term(x, ω, t)
    _fourier_term(x, ω, t, cos)
end


"""
$(TYPEDSIGNATURES)
Returns the coefficient of sin(ωt) in `x`.
"""
function fourier_sin_term(x, ω, t)
    _fourier_term(x, ω, t, sin)
end


function _fourier_term(x::Equation, ω, t, f)
    Equation(_fourier_term(x.lhs, ω, t, f) , _fourier_term(x.rhs, ω, t, f))
end

"Return the coefficient of f(ωt) in `x` where `f` is a cos or sin."
function _fourier_term(x, ω, t, f)
    term = x * f(ω * t)
    term = trig_reduce(term)
    indep = get_independent(term, t)
    ft = Num(simplify_complex(Symbolics.expand(indep)))
    ft = !isequal(ω, 0) ? 2*ft : ft # extra factor in case ω = 0 !
    Symbolics.expand(ft)
end

"Simplify fraction a/b + c/d = (ad + bc)/bd"
add_div(x) =  Num(Postwalk(add_with_div, similarterm=frac_similarterm)(unwrap(x)))

"Expand all sin/cos powers in `x`."
function trig_reduce(x)
    x = add_div(x) # a/b + c/d = (ad + bc)/bd
    x = expand(x) # open all brackets
    x = trig_to_exp(x)
    x = expand_all(x) # expand products of exponentials
    x = simplify_exp_products(x) # simplify products of exps
    x = exp_to_trig(x)
    x = Num(simplify_complex(expand(x)))
    simplify_fractions(x) # (a*c^2 + b*c)/c^2 = (a*c + b)/c
end


"Return the highest power of `y` occuring in the term `x`."
function max_power(x::Num, y::Num)
    terms = get_all_terms(x)
    powers = power_of.(terms, y)
    maximum(powers)
end

max_power(x::Vector{Num}, y::Num) = maximum(max_power.(x, y))
max_power(x::Complex, y::Num) = maximum(max_power.([x.re, x.im], y))
max_power(x, t) = max_power(Num(x), Num(t))

"Return the power of `y` in the term `x`"
function power_of(x::Num, y::Num)
    issym(y.val) ? nothing : error("power of " * string(y) * " is ambiguous")
    power_of(x.val, y.val)
end

function power_of(x::BasicSymbolic, y::BasicSymbolic)
    if ispow(x) && issym(y)
        return isequal(x.base, y) ? x.exp : 0
    elseif issym(x) && issym(y)
        return isequal(x, y) ? 1 : 0
    else
        return 0
    end
end

power_of(x, y) = 0
