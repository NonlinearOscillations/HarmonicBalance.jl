using Symbolics: unwrap, istree, operation, arguments, issym, diff2term
export ode_order_lowering, var_from_nested_derivative, lower_varname,first_order_transform!


function first_order_transform!(diff_eom::DifferentialEquation, time)
    eqs′, states′ = ode_order_lowering(diff_eom.equations, time, diff_eom.harmonics)
    diff_eom.equations = eqs′
    diff_eom.harmonics = states′
    nothing
end

function ode_order_lowering(equations, iv, harmonics)
    states = unwrap.(collect(keys(harmonics)))
    eqs = unwrap.(collect(values(equations)))

    var_order = OrderedDict{Any, Int}()
    D = Differential(iv)
    diff_eqs = similar(equations)
    diff_vars = similar(harmonics)

    for (i, eq) in enumerate(eqs)
        var, maxorder = var_from_nested_derivative(eq.lhs)
        maxorder > get(var_order, var, 1) && (var_order[var] = maxorder)

        var′ = lower_varname(var, iv, maxorder - 1)
        rhs′ = diff2term(eq.rhs)

        diff_vars[var′] = harmonics[var]
        diff_eqs[var′] =  D(var′) ~ rhs′
    end

    for (var, order) in var_order
        for o in (order - 1):-1:1
            lvar = lower_varname(var, iv, o - 1)
            rvar = lower_varname(var, iv, o)

            diff_vars[lvar] = harmonics[var]
            diff_eqs[lvar] = D(lvar) ~ rvar
        end
    end

    return (diff_eqs, diff_vars)
end

function var_from_nested_derivative(x,i=0)
    x = unwrap(x)
    if issym(x)
        (x, i)
    elseif istree(x)
        operation(x) isa Differential ?
            var_from_nested_derivative(first(arguments(x)), i + 1) : (x, i)
    else
        error("Not a well formed derivative expression $x")
    end
end

function lower_varname(var, idv, order)
    order == 0 && return var
    D = Differential(idv)
    for _ in 1:order
        var = D(var)
    end
    return diff2term(var)
end
