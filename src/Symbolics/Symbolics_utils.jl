"""
$(TYPEDSIGNATURES)

Perform substitutions in `rules` on `x`.
`include_derivatives=true` also includes all derivatives of the variables of the keys of `rules`.
"""
function substitute_all(
    x::T, rules::Dict; include_derivatives=true
)::T where {T<:Union{Equation,Num}}
    if include_derivatives
        rules = merge(
            rules,
            Dict([Differential(var) => Differential(rules[var]) for var in keys(rules)]),
        )
    end
    return substitute(x, rules)
end

"Variable substitution - dictionary"
function substitute_all(dict::Dict, rules::Dict)::Dict
    new_keys = substitute_all.(keys(dict), rules)
    new_values = substitute_all.(values(dict), rules)
    return Dict(zip(new_keys, new_values))
end

function substitute_all(
    v::Union{Array{Num},Array{Equation}}, rules::Union{Dict,Pair,Vector}
)
    return [substitute_all(x, rules) for x in v]
end
function substitute_all(x::Union{Num,Equation}, rules::Union{Pair,Vector,OrderedDict})
    return substitute_all(x, Dict(rules))
end
function substitute_all(x::Complex{Num}, rules::Union{Pair,Vector,OrderedDict,Dict})
    return substitute_all(Num(x.re.val.arguments[1]), rules)
end
substitute_all(x, rules) = substitute_all(Num(x), rules)


###
# STUFF BELOW IS MAINLY FOR FOURIER-TRANSFORMING
###

get_independent(x::Num, t::Num) = get_independent(x.val, t)
function get_independent(x::Complex{Num}, t::Num)
    return get_independent(x.re, t) + im * get_independent(x.im, t)
end
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
function get_all_terms(x::Equation)
    return unique(cat(get_all_terms(Num(x.lhs)), get_all_terms(Num(x.rhs)); dims=1))
end

_get_all_terms_mul(x) = Num.(SymbolicUtils.arguments(x))
_get_all_terms_div(x) = Num.([_get_all_terms(x.num)..., _get_all_terms(x.den)...])
_get_all_terms(x) = Num(x)

function _get_all_terms_add(x)::Vector{Num}
    list = []
    for term in keys(x.dict)
        list = cat(list, _get_all_terms(term); dims=1)
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


"Return true if `f` is a function of `var`."
is_function(f, var) = any(isequal.(get_variables(f), var))


"Simplify fraction a/b + c/d = (ad + bc)/bd"
add_div(x) = Num(Postwalk(add_with_div; maketerm=frac_maketerm)(unwrap(x)))

"Return the highest power of `y` occuring in the term `x`."
function max_power(x::Num, y::Num)
    terms = get_all_terms(x)
    powers = power_of.(terms, y)
    return maximum(powers)
end

max_power(x::Vector{Num}, y::Num) = maximum(max_power.(x, y))
max_power(x::Complex, y::Num) = maximum(max_power.([x.re, x.im], y))
max_power(x, t) = max_power(Num(x), Num(t))

"Return the power of `y` in the term `x`"
function power_of(x::Num, y::Num)
    issym(y.val) ? nothing : error("power of " * string(y) * " is ambiguous")
    return power_of(x.val, y.val)
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
