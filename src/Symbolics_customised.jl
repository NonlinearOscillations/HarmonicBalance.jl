import Symbolics.SymbolicUtils: quick_cancel; export quick_cancel
import Symbolics.SymbolicUtils: Postwalk


# change SymbolicUtils' quick_cancel to simplify powers of fractions correctly
function quick_cancel(x::Term, y::Term)
	if x.f == exp && y.f == exp
		return exp(x.arguments[1] - y.arguments[1]), 1
	else
		return x,y
	end
end

quick_cancel(x::Term, y::Pow) = y.base isa Term && y.base.f == exp ? quick_cancel(x, expand_exp_power(y)) : x,y


"Returns true if expr is an exponential"
is_exp(expr) = isterm(expr) && expr.f == exp

"Expand powers of exponential such that exp(x)^n => exp(x*n) "
expand_exp_power(expr) = ispow(expr) && is_exp(expr.base) ? exp(expr.base.arguments[1] * expr.exp) : expr
expand_exp_power_add(expr) = sum([expand_exp_power(arg) for arg in arguments(expr)])
expand_exp_power_mul(expr) = prod([expand_exp_power(arg) for arg in arguments(expr)])
expand_exp_power(expr::Num) = expand_exp_power(expr.val)

function expand_exp_power(expr::BasicSymbolic)
    if isadd(expr)
        return expand_exp_power_add(expr)
    elseif ismul(expr)
        return expand_exp_power_mul(expr)
    else
        return ispow(expr) && is_exp(expr.base) ? exp(expr.base.arguments[1] * expr.exp) : expr
    end
end

"Expands using SymbolicUtils.expand and expand_exp_power (changes exp(x)^n to exp(x*n)"
expand_all(x) = Postwalk(expand_exp_power)(SymbolicUtils.expand(x))
expand_all(x::Complex{Num}) = expand_all(x.re) + im* expand_all(x.im)
expand_all(x::Num) = Num(expand_all(x.val))


"Apply a function f on every member of a sum or a product"
_apply_termwise(f, x) = f(x)
function _apply_termwise(f, x::BasicSymbolic)
    if isadd(x)
        return sum([f(arg) for arg in arguments(x)])
    elseif ismul(x)
        return prod([f(arg) for arg in arguments(x)])
    elseif isdiv(x)
        return _apply_termwise(f, x.num) / _apply_termwise(f, x.den)
    else
        return  f(x)
    end
end

simplify_complex(x::Complex) = isequal(x.im, 0) ? x.re : x.re + im*x.im
simplify_complex(x) = x
function simplify_complex(x::BasicSymbolic)
    if isadd(x) || ismul(x) ||  isdiv(x)
        return _apply_termwise(simplify_complex, x)
    else
        return x
    end
end


"Simplify products of exponentials such that exp(a)*exp(b) => exp(a+b)
This is included in SymbolicUtils as of 17.0 but the method here avoid other simplify calls"
function simplify_exp_products_mul(expr)
    ind = findall(x -> is_exp(x), arguments(expr))
    rest_ind = setdiff(1:length(arguments(expr)), ind)
    rest = isempty(rest_ind) ? 1 : prod(arguments(expr)[rest_ind])
    total = isempty(ind) ? 0 : sum(getindex.(arguments.(arguments(expr)[ind]), 1))
    SymbolicUtils.is_literal_number(total) ? (total == 0 && return rest) : return rest * exp(total)
end

simplify_exp_products(x::Complex{Num}) = Complex{Num}(simplify_exp_products(x.re.val), simplify_exp_products(x.im.val))
simplify_exp_products(x::Num) = simplify_exp_products(x.val)
simplify_exp_products(x) = x

function simplify_exp_products(expr::BasicSymbolic)
    if isadd(expr) || isdiv(expr)
        return _apply_termwise(simplify_exp_products, expr)
    elseif ismul(expr)
        return simplify_exp_products_mul(expr)
    else
        return expr
    end
end

function exp_to_trig(x::BasicSymbolic)
    if isadd(x) || isdiv(x) || ismul(x)
        return _apply_termwise(exp_to_trig, x)
    elseif isterm(x) && x.f == exp
        arg = first(x.arguments)
        trigarg = Symbolics.expand(-im*arg) # the argument of the to-be trig function
        trigarg = simplify_complex(trigarg)

        # put arguments of trigs into a standard form such that sin(x) = -sin(-x), cos(x) = cos(-x) are recognized
        if isadd(trigarg)
            first_symbol = minimum(cat(string.(arguments(trigarg)), string.(arguments(-trigarg)), dims=1))

            # put trigarg => -trigarg the lowest alphabetic argument of trigarg is lower than that of -trigarg
            # this is a meaningless key but gives unique signs to all sums
            is_first = minimum(string.(arguments(trigarg))) == first_symbol
            return is_first ? cos(-trigarg) -im*sin(-trigarg) : cos(trigarg)+im* sin(trigarg)
        end
        return ismul(trigarg) && trigarg.coeff < 0 ? cos(-trigarg) -im*sin(-trigarg) : cos(trigarg)+im* sin(trigarg)
    else
        return x
    end
end

exp_to_trig(x)=x
exp_to_trig(x::Num) = exp_to_trig(x.val)
exp_to_trig(x::Complex{Num}) = exp_to_trig(x.re)+im* exp_to_trig(x.im)


# sometimes, expressions get stored as Complex{Num} with no way to decode what real(x) and imag(x)
# this overloads the Num constructor to return a Num if x.re and x.im have similar arguments
Num(x::Complex{Num})::Num = isequal(x.re.val.arguments, x.im.val.arguments) ? Num(first(x.re.val.arguments)) : error("Cannot convert Complex{Num} " * string(x) * " to Num")


#=
chop(x) = x
chop(x::Complex{Int64})= Int64(x)
chop(x::Complex{Float64}) = Float64(x)
chop(x::Add) = _apply_termwise(chop, x)
chop(x::Mul) = _apply_termwise(chop, x)
chop(x::Num) = chop(x.val)
chop(x::Complex{Num}) = Complex{Num}(x.re, x.im)


#expand_fraction(expr::Div) = expr.num isa Add ? sum([arg / expr.den for arg in arguments(expr.num)]) : expr
#expand_fraction(expr) = expr

#simplify_parts(expr::Num) = simplify_parts(expr.val)
#simplify_parts(expr::Add) = sum([simplify_fractions(arg) for arg in arguments(expr)])
#simplify_parts(expr) = simplify_fractions(expr)
=#
