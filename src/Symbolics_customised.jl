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
    add_with_div,
    frac_maketerm,
    @rule,
    @acrule,
    @compactified,
    isliteral,
    Chain
using SymbolicUtils.TermInterface: issym
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

# ∨ Complex{Num} can contain
is_false_or_zero(x) = x === false || _iszero(x)
is_literal_complex(x) = isliteral(Complex)(x)
is_literal_real(x) = isliteral(Real)(x)
function is_not_complex(x)
    return !is_literal_real(x) && (is_literal_complex(x) && is_false_or_zero(unwrap(x.im)))
end

"Expands using SymbolicUtils.expand and expand_exp_power (changes exp(x)^n to exp(x*n)"
const expand_exp_power = @rule(exp(~x)^(~y) => exp(~x * ~y))
expand_all(x) = simplify(expand(x); rewriter=Postwalk(expand_exp_power))
expand_all(x::Num) = wrap(expand_all(unwrap(x)))
function expand_all(x::Complex{Num})
    re_val = is_false_or_zero(unwrap(x.re)) ? 0.0 : expand_all(x.re)
    im_val = is_false_or_zero(unwrap(x.im)) ? 0.0 : expand_all(x.im)
    return re_val + im * im_val
end # This code is stupid, we can just use simplify

simplify_complex(x) = Postwalk(@rule(~x::is_not_complex => real(~x)))(x)

"Simplify products of exponentials such that exp(a)*exp(b) => exp(a+b)"
const simplify_exp_mul = @acrule(exp(~x) * exp(~y) => _iszero(~x + ~y) ? 1 : exp(~x + ~y))
simplify_exp_products(x) = Postwalk(simplify_exp_mul)(x)
function simplify_exp_products(x::Complex{Num})
    re_val = is_false_or_zero(unwrap(x.re)) ? 0.0 : simplify_exp_products(x.re)
    im_val = is_false_or_zero(unwrap(x.im)) ? 0.0 : simplify_exp_products(x.im)
    return re_val + im * im_val
end

const sin_euler = @rule(sin(~x) => (exp(im * ~x) - exp(-im * ~x)) / (2 * im))
const cos_euler = @rule(cos(~x) => (exp(im * ~x) + exp(-im * ~x)) / 2)
trig_to_exp(x) = simplify(x; rewriter=Postwalk(Chain([sin_euler, cos_euler])))
trig_to_exp(x::Num) = wrap(trig_to_exp(unwrap(x)))
function trig_to_exp(x::Complex{Num})
    re_val = is_false_or_zero(unwrap(x.re)) ? 0.0 : trig_to_exp(x.re)
    im_val = is_false_or_zero(unwrap(x.im)) ? 0.0 : trig_to_exp(x.im)
    return re_val + im * im_val
end

# @compactified is what SymbolicUtils uses internally
# function _apply_termwise(f, x::BasicSymbolic)
#     @compactified x::BasicSymbolic begin
#     Add  => sum([f(arg) for arg in arguments(x)])
#     Mul  => prod([f(arg) for arg in arguments(x)])
#     Div  =>  _apply_termwise(f, x.num) / _apply_termwise(f, x.den)
#     _    => f(x)
#     end
# end

# TODO ∨ try symbolicUtils.substitute such that maketerm gets called
reparse(x) = parse_expr_to_symbolic(Meta.parse(string(x)), @__MODULE__)
# ^ parsing and reevaluting makes it that (1 - 0.0im) becomes (1)

function make_positive_trig(x::Num)
    all_terms = get_all_terms(x)
    trigs = filter(z -> is_trig(z), all_terms)

    rules = []
    for trig in trigs
        is_pow = ispow(trig.val) # trig is either a trig or a power of a trig
        power = is_pow ? trig.val.exp : 1
        arg = is_pow ? arguments(trig.val.base)[1] : arguments(trig.val)[1]
        type = is_pow ? operation(trig.val.base) : operation(trig.val)
        negative =
            !issym(arg) && prod(Number.(filter(x -> x isa Number, arguments(arg)))) < 0

        if negative
            if type == cos
                term = cos(-arg)
            elseif type == sin
                term = (-1)^power * sin(-arg)^power
            end
            append!(rules, [trig => term])
        end
    end
    result = Symbolics.substitute(x, Dict(rules))
    return result
end

function exp_to_trig(z::Complex{Num})
    z = reparse(z)
    return simplify_complex(make_positive_trig(z.re) + im * make_positive_trig(z.im))
end
