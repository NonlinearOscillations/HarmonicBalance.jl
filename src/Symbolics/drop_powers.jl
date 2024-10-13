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
    Symbolics.@variables ϵ
    subs_expr = deepcopy(expr)
    rules = Dict([var => ϵ * var for var in unique(vars)])
    subs_expr = Symbolics.expand(substitute_all(subs_expr, rules))
    max_deg = max_power(subs_expr, ϵ)
    removal = Dict([ϵ^d => Num(0) for d in deg:max_deg])
    res = substitute_all(substitute_all(subs_expr, removal), Dict(ϵ => Num(1)))
    return Symbolics.expand(res)
    #res isa Complex ? Num(res.re.val.arguments[1]) : res
end

function drop_powers(expr::Vector{Num}, var::Vector{Num}, deg::Int)
    return [drop_powers(x, var, deg) for x in expr]
end

# calls the above for various types of the first argument
function drop_powers(eq::Equation, var::Vector{Num}, deg::Int)
    return drop_powers(eq.lhs, var, deg) .~ drop_powers(eq.lhs, var, deg)
end
function drop_powers(eqs::Vector{Equation}, var::Vector{Num}, deg::Int)
    return [
        Equation(drop_powers(eq.lhs, var, deg), drop_powers(eq.rhs, var, deg)) for eq in eqs
    ]
end
drop_powers(expr, var::Num, deg::Int) = drop_powers(expr, [var], deg)
drop_powers(x, vars, deg::Int) = drop_powers(Num(x), vars, deg)
