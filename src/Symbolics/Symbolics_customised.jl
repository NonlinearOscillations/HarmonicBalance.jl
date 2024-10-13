
expand_all(x::Num) = Num(expand_all(x.val))
_apply_termwise(f, x::Num) = wrap(_apply_termwise(f, unwrap(x)))

"Expands using SymbolicUtils.expand and expand_exp_power (changes exp(x)^n to exp(x*n)"
function expand_all(x)
    result = Postwalk(expand_exp_power)(SymbolicUtils.expand(x))
    return isnothing(result) ? x : result
end
expand_all(x::Complex{Num}) = expand_all(x.re) + im * expand_all(x.im)

"Apply a function f on every member of a sum or a product"
function _apply_termwise(f, x::BasicSymbolic)
    @compactified x::BasicSymbolic begin
    Add  => sum([f(arg) for arg in arguments(x)])
    Mul  => prod([f(arg) for arg in arguments(x)])
    Div  =>  _apply_termwise(f, x.num) / _apply_termwise(f, x.den)
    _    => f(x)
    end
end

simplify_complex(x::Complex) = isequal(x.im, 0) ? x.re : x.re + im * x.im
simplify_complex(x) = x
function simplify_complex(x::BasicSymbolic)
    if isadd(x) || ismul(x) || isdiv(x)
        return _apply_termwise(simplify_complex, x)
    else
        return x
    end
end



# sometimes, expressions get stored as Complex{Num} with no way to decode what real(x) and imag(x)
# this overloads the Num constructor to return a Num if x.re and x.im have similar arguments
function Num(x::Complex{Num})::Num
    if x.re.val isa Float64 && x.im.val isa Float64
        return Num(x.re.val)
    else
        if isequal(x.re.val.arguments, x.im.val.arguments)
            Num(first(x.re.val.arguments))
        else
            error("Cannot convert Complex{Num} " * string(x) * " to Num")
        end
    end
end
# ^ This function commits type-piracy with Symbolics.jl. We should change this.
