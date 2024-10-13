using Symbolics
using SymbolicUtils:
    SymbolicUtils,
    Postwalk,
    Sym,
    BasicSymbolic,
    isterm,
    ispow,
    isadd,
    isdiv,
    ismul,
    issym,
    add_with_div,
    frac_maketerm,
    @rule,
    @acrule,
    @compactified,
    isliteral,
    Chain
using Symbolics:
    Symbolics,
    Num,
    unwrap,
    wrap,
    get_variables,
    simplify,
    expand_derivatives,
    build_function,
    Equation,
    Differential,
    @variables,
    arguments,
    simplify_fractions,
    substitute,
    term,
    expand,
    operation,
    _iszero

const expand_exp_power = @rule(exp(~x)^(~y) => exp(~x * ~y))
const simplify_exp_mul = @acrule(exp(~x) * exp(~y) => _iszero(~x + ~y) ? 1 : exp(~x + ~y))
const sin_euler = @rule(sin(~x) => (exp(im*~x)-exp(-im*~x))/(2*im))
const cos_euler = @rule(cos(~x) => (exp(im*~x)+exp(-im*~x))/2)

# ∨ Complex{Num} can contain
is_false_or_zero(x) = x === false || _iszero(x)
is_literal_complex(x) = isliteral(Complex)(x)
is_literal_real(x) = isliteral(Real)(x)
function is_not_complex(x)
    return !is_literal_real(x) && (is_literal_complex(x) && is_false_or_zero(unwrap(x.im)))
end

"Expands using SymbolicUtils.expand and expand_exp_power (changes exp(x)^n to exp(x*n)"

expand_all(x) = simplify(expand(x), rewriter=Postwalk(expand_exp_power))
expand_all(x::Num) = wrap(expand_all(unwrap(x)))
function expand_all(x::Complex{Num})
    re_val = is_false_or_zero(unwrap(x.re)) ? 0.0 : expand_all(x.re)
    im_val = is_false_or_zero(unwrap(x.im)) ? 0.0 : expand_all(x.im)
    return re_val + im * im_val
end # This code is stupid, we can just use simplify

simplify_complex(x) = Postwalk(@rule(~x::is_not_complex => real(~x)))(x)

"Simplify products of exponentials such that exp(a)*exp(b) => exp(a+b)"
simplify_exp_products(x) = Postwalk(simplify_exp_mul)(x)
function simplify_exp_products(x::Complex{Num})
    re_val = is_false_or_zero(unwrap(x.re)) ? 0.0 : simplify_exp_products(x.re)
    im_val = is_false_or_zero(unwrap(x.im)) ? 0.0 : simplify_exp_products(x.im)
    return re_val + im * im_val
end

trig_to_exp(x) = simplify(x; rewriter= Postwalk(Chain([sin_euler,cos_euler])))
trig_to_exp(x::Num) = wrap(trig_to_exp(unwrap(x)))
function trig_to_exp(x::Complex{Num})
    re_val = is_false_or_zero(unwrap(x.re)) ? 0.0 : trig_to_exp(x.re)
    im_val = is_false_or_zero(unwrap(x.im)) ? 0.0 : trig_to_exp(x.im)
    return re_val + im * im_val
end # This code is stupid, we can just use simplify


# function exp_to_trig(x::BasicSymbolic)
#     if isadd(x) || isdiv(x) || ismul(x)
#         return _apply_termwise(exp_to_trig, x)
#     elseif isterm(x) && x.f == exp
#         arg = first(x.arguments)
#         trigarg = simplify_complex(Symbolics.expand(-im * arg))

#         return isadd(trigarg) ? handle_add_trigarg(trigarg) : handle_prod_trigarg(trigarg)
#     else
#         return x
#     end
# end
# exp_to_trig(x) = x
# exp_to_trig(x::Num) = Num(exp_to_trig(x.val))
# exp_to_trig(x::Complex{Num}) = exp_to_trig(x.re) + im * exp_to_trig(x.im)

# is_literal_complex(x) = isliteral(Complex)(x)
# euler = @rule(exp(im*~~z) => exp(real(prod(~~z)))*(cos(imag(prod(~~z))) + im * sin(imag(prod(~~z)))))
# euler_identity(x) = Postwalk(euler)(x)
# euler(exp(im*a))
# euler_identity(exp(im*a))
# simplify_complex(z + a) |> typeof

# euler(arg, f=+) = f(cos(arg), im * sin(arg))
# function handle_add_trigarg(trigarg)
#     terms = string.(arguments(trigarg))
#     min_terms = string.(arguments(-trigarg))
#     symbols = cat(terms, min_terms; dims=1)

#     first_symbol = minimum(symbols)
#     is_first = minimum(string.(arguments(trigarg))) == first_symbol
#     return is_first ? euler(-trigarg, -) : euler(trigarg)
# end

# function handle_non_add_trigarg(trigarg)
#     return (ismul(trigarg) && trigarg.coeff < 0) ? euler(-trigarg, -) : euler(trigarg)
# end

# sometimes, expressions get stored as Complex{Num} with no way to decode what real(x) and imag(x)
# this overloads the Num constructor to return a Num if x.re and x.im have similar arguments
# function Num(x::Complex{Num})::Num
#     if x.re.val isa Float64 && unwrap(x.im.val) isa Float64
#         return Num(x.re.val)
#     else
#         if isequal(x.re.val.arguments, unwrap(x.im.val).arguments)
#             Num(first(x.re.val.arguments))
#         else
#             error("Cannot convert Complex{Num} " * string(x) * " to Num")
#         end
#     end
# end
# ^ This function commits type-piracy with Symbolics.jl. We should change this.
