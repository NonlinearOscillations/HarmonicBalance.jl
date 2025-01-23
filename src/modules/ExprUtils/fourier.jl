
"Expand all sin/cos powers in `x`."
function trig_reduce(x)
    x = add_div(x) # a/b + c/d = (ad + bc)/bd
    x = expand(x) # open all brackets
    x = trig_to_exp(x)
    x = expand_all(x) # expand products of exponentials
    x = simplify_exp_products(x) # simplify products of exps
    x = exp_to_trig(x)
    x = Num(simplify_complex(expand(x)))
    return x# simplify_fractions(x)# (a*c^2 + b*c)/c^2 = (a*c + b)/c
end

"Return true if `f` is a sin or cos."
function is_trig(f::Num)
    f = ispow(f.val) ? f.val.base : f.val
    isterm(f) && SymbolicUtils.operation(f) ∈ [cos, sin] && return true
    return false
end

"""
$(TYPEDSIGNATURES)
Returns the coefficient of cos(ωt) in `x`.
"""
function fourier_cos_term(x, ω, t)
    return _fourier_term(x, ω, t, cos)
end

"Simplify fraction a/b + c/d = (ad + bc)/bd"
add_div(x) = wrap(Postwalk(add_with_div; maketerm=frac_maketerm)(unwrap(x)))

"""
$(TYPEDSIGNATURES)
Returns the coefficient of sin(ωt) in `x`.
"""
function fourier_sin_term(x, ω, t)
    return _fourier_term(x, ω, t, sin)
end

function _fourier_term(x::Equation, ω, t, f)
    return Equation(_fourier_term(x.lhs, ω, t, f), _fourier_term(x.rhs, ω, t, f))
end

"Return the coefficient of f(ωt) in `x` where `f` is a cos or sin."
function _fourier_term(x, ω, t, f)
    term = x * f(ω * t)
    term = trig_reduce(term)
    indep = get_independent(term, t)
    ft = Num(simplify_complex(Symbolics.expand(indep)))
    ft = !isequal(ω, 0) ? 2 * ft : ft # extra factor in case ω = 0 !
    return Symbolics.expand(ft)
end

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
            term = Complex{Num}((exp(im * arg) + exp(-im * arg))^power * (1//2)^power, 0)
        elseif type == sin
            term =
                (1 * im^power) *
                Complex{Num}(((exp(-im * arg) - exp(im * arg)))^power * (1//2)^power, 0)
        end
        # avoid Complex{Num} where possible as this causes bugs
        # instead, the Nums store SymbolicUtils Complex types
        term = Num(Symbolics.expand(term.re.val + im * term.im.val))
        append!(rules, [trig => term])
    end

    result = Symbolics.substitute(x, Dict(rules))
    return convert_to_Num(result)
end
convert_to_Num(x::Complex{Num})::Num = Num(first(x.re.val.arguments))
convert_to_Num(x::Num)::Num = x

function exp_to_trig(x::BasicSymbolic)
    if isadd(x) || isdiv(x) || ismul(x)
        return _apply_termwise(exp_to_trig, x)
    elseif isterm(x) && x.f == exp
        arg = first(x.arguments)
        trigarg = Symbolics.expand(-im * arg) # the argument of the to-be trig function
        trigarg = simplify_complex(trigarg)

        # put arguments of trigs into a standard form such that sin(x) = -sin(-x), cos(x) = cos(-x) are recognized
        if isadd(trigarg)
            first_symbol = minimum(
                cat(string.(arguments(trigarg)), string.(arguments(-trigarg)); dims=1)
            )

            # put trigarg => -trigarg the lowest alphabetic argument of trigarg is lower than that of -trigarg
            # this is a meaningless key but gives unique signs to all sums
            is_first = minimum(string.(arguments(trigarg))) == first_symbol
            return if is_first
                cos(-trigarg) - im * sin(-trigarg)
            else
                cos(trigarg) + im * sin(trigarg)
            end
        end
        return if ismul(trigarg) && trigarg.coeff < 0
            cos(-trigarg) - im * sin(-trigarg)
        else
            cos(trigarg) + im * sin(trigarg)
        end
    else
        return x
    end
end

exp_to_trig(x) = x
exp_to_trig(x::Num) = exp_to_trig(x.val)
exp_to_trig(x::Complex{Num}) = exp_to_trig(x.re) + im * exp_to_trig(x.im)
