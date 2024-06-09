export first_order_transform!, is_rearranged_standard, rearrange_standard!, get_equations

get_equations(eom::DifferentialEquation) = collect(values(eom.equations))

# TODO: check the degree of the eom
function is_rearranged_standard(eom::DifferentialEquation, degree = 2)
    tvar = get_independent_variables(eom)[1]
    D = Differential(tvar)^degree
    isequal(getfield.(values(eom.equations), :lhs), D.(get_variables(eom)))
end

function rearrange_standard!(eom::DifferentialEquation, degree = 2)
    tvar = get_independent_variables(eom)[1]
    D = Differential(tvar)^degree
    dvars = D.(get_variables(eom))
    rearrange!(eom, dvars)
end

function HarmonicBalance.rearrange!(eom::DifferentialEquation, new_lhs::Vector{Num})
    soln = Symbolics.solve_for(get_equations(eom), new_lhs, simplify = false, check = true)
    eom.equations = OrderedDict(zip(get_variables(new_lhs), new_lhs .~ soln))
    return nothing
end

function HarmonicBalance.rearrange(eom::DifferentialEquation, new_lhs::Vector{Num})
    new_eom = deepcopy(eom)
    rearrange!(new_eom, new_lhs)
    return new_eom
end

function first_order_transform!(diff_eom::DifferentialEquation, time)
    eqs′, states′ = ode_order_lowering(diff_eom.equations, time, diff_eom.harmonics)
    diff_eom.equations = eqs′
    diff_eom.harmonics = states′
    nothing
end

# taken from ModelingToolkit.jl
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
        diff_eqs[var′] = D(var′) ~ rhs′
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
