
expand_all(x::Num) = Num(expand_all(x.val))
_apply_termwise(f, x::Num) = wrap(_apply_termwise(f, unwrap(x)))

"Expands using SymbolicUtils.expand and expand_exp_power (changes exp(x)^n to exp(x*n)"
function expand_all(x)
    result = Postwalk(expand_exp_power)(SymbolicUtils.expand(x))
    return isnothing(result) ? x : result
end
expand_all(x::Complex{Num}) = expand_all(x.re) + im * expand_all(x.im)

function expand_fraction(x::BasicSymbolic)
    @compactified x::BasicSymbolic begin
        Add => _apply_termwise(expand_fraction, x)
        Mul => _apply_termwise(expand_fraction, x)
        Div => sum([arg / x.den for arg in arguments(x.num)])
        _   => x
    end
end
expand_fraction(x::Num) = Num(expand_fraction(x.val))

"Apply a function f on every member of a sum or a product"
function _apply_termwise(f, x::BasicSymbolic)
    @compactified x::BasicSymbolic begin
        Add => sum([f(arg) for arg in arguments(x)])
        Mul => prod([f(arg) for arg in arguments(x)])
        Div => _apply_termwise(f, x.num) / _apply_termwise(f, x.den)
        _   => f(x)
    end
end

simplify_complex(x::Complex) = isequal(x.im, 0) ? x.re : x.re + im * x.im
simplify_complex(x) = x
function simplify_complex(x::BasicSymbolic)
    @compactified x::BasicSymbolic begin
        Add => _apply_termwise(simplify_complex, x)
        Mul => _apply_termwise(simplify_complex, x)
        Div => _apply_termwise(simplify_complex, x)
        _   => x
    end
end

"""
$(TYPEDSIGNATURES)

Perform substitutions in `rules` on `x`.
`include_derivatives=true` also includes all derivatives of the variables of the keys of `rules`.
"""
Subtype = Union{Num,Equation,BasicSymbolic}
function substitute_all(x::Subtype, rules::Dict; include_derivatives=true)
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
Collections = Union{Dict,Pair,Vector,OrderedDict}
substitute_all(v::AbstractArray, rules) = [substitute_all(x, rules) for x in v]
substitute_all(x::Subtype, rules::Collections) = substitute_all(x, Dict(rules))
function substitute_all(x::Complex{Num}, rules::Collections)
    return substitute_all(x.re, rules) + im * substitute_all(x.im, rules)
end

get_independent(x::Num, t::Num) = get_independent(x.val, t)
function get_independent(x::Complex{Num}, t::Num)
    return get_independent(x.re, t) + im * get_independent(x.im, t)
end
get_independent(v::Vector{Num}, t::Num) = [get_independent(el, t) for el in v]
get_independent(x, t::Num) = x

function get_independent(x::BasicSymbolic, t::Num)
    @compactified x::BasicSymbolic begin
        Add  => sum([get_independent(arg, t) for arg in arguments(x)])
        Mul  => prod([get_independent(arg, t) for arg in arguments(x)])
        Div  => !is_function(x.den, t) ? get_independent(x.num, t) / x.den : 0
        Pow  => !is_function(x.base, t) && !is_function(x.exp, t) ? x : 0
        Term => !is_function(x, t) ? x : 0
        Sym  => !is_function(x, t) ? x : 0
        _    => x
    end
end

"Return all the terms contained in `x`"
get_all_terms(x::Num) = unique(_get_all_terms(Symbolics.expand(x).val))
function get_all_terms(x::Equation)
    return unique(cat(get_all_terms(Num(x.lhs)), get_all_terms(Num(x.rhs)); dims=1))
end
function _get_all_terms(x::BasicSymbolic)
    @compactified x::BasicSymbolic begin
        Add => vcat([_get_all_terms(term) for term in SymbolicUtils.arguments(x)]...)
        Mul => Num.(SymbolicUtils.arguments(x))
        Div => Num.([_get_all_terms(x.num)..., _get_all_terms(x.den)...])
        _   => Num(x)
    end
end
_get_all_terms(x) = Num(x)

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

"""
Counts the number of derivatives of a symbolic variable.
"""
function count_derivatives(x::Symbolics.BasicSymbolic)
    (Symbolics.isterm(x) || Symbolics.issym(x)) ||
        error("The input is not a single term or symbol")
    bool = Symbolics.is_derivative(x)
    return bool ? 1 + count_derivatives(first(arguments(x))) : 0
end
count_derivatives(x::Num) = count_derivatives(Symbolics.unwrap(x))
