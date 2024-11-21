"""
$(TYPEDEF)

Holds differential equation(s) of motion and a set of harmonics to expand each variable.
This is the primary input for `HarmonicBalance.jl`. After inputting the equations,
the harmonics ansatz needs to be specified using `add_harmonic!`.

# Fields
$(TYPEDFIELDS)

## Example
```julia-repl
julia> @variables t, x(t), y(t), ω0, ω, F, k;

# equivalent ways to enter the simple harmonic oscillator
julia> DifferentialEquation(d(x,t,2) + ω0^2 * x - F * cos(ω*t), x);
julia> DifferentialEquation(d(x,t,2) + ω0^2 * x ~ F * cos(ω*t), x);

# two coupled oscillators, one of them driven
julia> DifferentialEquation(
    [d(x,t,2) + ω0^2 * x - k*y, d(y,t,2) + ω0^2 * y - k*x] .~ [F * cos(ω*t), 0], [x,y]
);
```
"""
mutable struct DifferentialEquation
    """Assigns to each variable an equation of motion."""
    equations::OrderedDict{Num,Equation}
    """Assigns to each variable a set of harmonics."""
    harmonics::OrderedDict{Num,OrderedSet{Num}}

    function DifferentialEquation(eqs)
        return new(eqs, OrderedDict(var => OrderedSet() for var in keys(eqs)))
    end

    # uses the above constructor if no harmonics defined
    function DifferentialEquation(eqs::Vector{Equation}, vars::Vector{Num})
        return DifferentialEquation(OrderedDict(zip(vars, eqs)))
    end

    # if expressions are entered instead of equations, automatically set them = 0
    function DifferentialEquation(exprs::Vector{Num}, vars::Vector{Num})
        return DifferentialEquation(exprs .~ Int(0), vars)
    end

    function DifferentialEquation(eq::Equation, var::Num)
        typerhs = typeof(eq.rhs)
        typelhs = typeof(eq.lhs)
        if eq.rhs isa AbstractVector || eq.lhs isa AbstractVector
            throw(
                ArgumentError(
                    "The equation is of the form $(typerhs)~$(typelhs) is not supported.
                    Commenly one forgot to broadcast the equation symbol `~`."
                ),
            )
        end
        return DifferentialEquation([eq], [var])
    end
    function DifferentialEquation(eq::Equation, vars::Vector{Num})
        typerhs = typeof(eq.rhs)
        typelhs = typeof(eq.lhs)
        throw(
            ArgumentError(
                "The variables are of type $(typeof(vars)). Whereas you gave one equation of
                type $(typerhs)~$(typelhs). Commenly one forgot to broadcast the equation symbol `~`.",
            ),
        )
    end
    DifferentialEquation(lhs::Num, var::Num) = DifferentialEquation([lhs ~ Int(0)], [var])
end

"show method of the type `DifferentialEquation`"
function Base.show(io::IO, diff_eq::DifferentialEquation)
    println(io, "System of ", length(keys(diff_eq.equations)), " differential equations")
    println(io, "Variables:       ", join(keys(diff_eq.equations), ", "))
    print(io, "Harmonic ansatz: ")
    for var in keys(diff_eq.harmonics)
        print(io, string(var), " => ", join(string.(diff_eq.harmonics[var]), ", "))
        print(io, ";   ")
    end
    println(io, "\n")
    return [println(io, eq) for eq in values(diff_eq.equations)]
end
"
Displays the fields of the differential equation object.
"
Base.show(eom::DifferentialEquation) = show_fields(eom)

"""
$(TYPEDSIGNATURES)

Add the harmonic `ω` to the harmonic ansatz used to expand the variable `var` in `diff_eom`.

## Example

# define the simple harmonic oscillator and specify that x(t) oscillates with frequency ω
```julia-repl
julia> @variables t, x(t), y(t), ω0, ω, F, k;
julia> diff_eq = DifferentialEquation(d(x,t,2) + ω0^2 * x ~ F * cos(ω*t), x);
julia> add_harmonic!(diff_eq, x, ��) # expand x using ω

System of 1 differential equations
Variables:       x(t)
Harmonic ansatz: x(t) => ω;

(ω0^2)*x(t) + Differential(t)(Differential(t)(x(t))) ~ F*cos(t*ω)
```
"""
function add_harmonic!(diff_eom::DifferentialEquation, var::Num, ω)
    push!.(Ref(diff_eom.harmonics[var]), ω)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the dependent variables of `diff_eom`.
"""
Symbolics.get_variables(diff_eom::DifferentialEquation)::Vector{Num} =
    collect(keys(diff_eom.equations))

"""
$(TYPEDSIGNATURES)

Check if all equations in `diff_eom` are harmonic with respect to `t`. The function takes a
differential equation system `diff_eom` and a variable `t`, and returns `true` if all equations
are harmonic with respect to `t`, otherwise it returns `false`.
"""
ExprUtils.is_harmonic(diff_eom::DifferentialEquation, t::Num)::Bool =
    all([is_harmonic(eq, t) for eq in values(diff_eom.equations)])

"""
$(TYPEDSIGNATURES)

Return the independent dependent variables of `diff_eom`.
"""
function get_independent_variables(diff_eom::DifferentialEquation)
    return Num.(flatten(unique([x.val.arguments for x in keys(diff_eom.equations)])))
end


get_equations(eom::DifferentialEquation) = collect(values(eom.equations))

"""
$(TYPEDSIGNATURES)

Checks if the differential equations in `eom` are arranged in standard form, where the highest
derivative of each variable appears isolated on the left-hand side. The default degree is 2,
corresponding to second-order differential equations.
"""
function is_rearranged_standard(eom::DifferentialEquation, degree=2)
    tvar = get_independent_variables(eom)[1]
    D = Differential(tvar)^degree
    return isequal(getfield.(values(eom.equations), :lhs), D.(get_variables(eom)))
end

"""
$(TYPEDSIGNATURES)

Rearranges the differential equations in `eom` to standard form, where the highest derivative
of each variable (specified by `degree`, default 2) appears isolated on the left-hand side.
Modifies the equations in place.
"""
function rearrange_standard!(eom::DifferentialEquation, degree=2)
    tvar = get_independent_variables(eom)[1]
    D = Differential(tvar)^degree
    dvars = D.(get_variables(eom))
    return rearrange!(eom, dvars)
end

"""
$(TYPEDSIGNATURES)

Rearranges the equations in `eom` such that the expressions in `new_lhs` appear isolated on
the left-hand sides. Uses symbolic linear solving to determine the right-hand sides. Modifies
the equations in place.
"""
function rearrange!(eom::DifferentialEquation, new_lhs::Vector{Num})
    soln = Symbolics.symbolic_linear_solve(
        get_equations(eom), new_lhs; simplify=false, check=true
    )
    eom.equations = OrderedDict(zip(get_variables_nums(new_lhs), new_lhs .~ soln))
    return nothing
end

"""
$(TYPEDSIGNATURES)

Creates a new differential equation system by rearranging the equations in `eom` such that
the expressions in `new_lhs` appear isolated on the left-hand sides. Similar to `rearrange!`
but returns a new system instead of modifying in place.
"""
function rearrange(eom::DifferentialEquation, new_lhs::Vector{Num})
    new_eom = deepcopy(eom)
    rearrange!(new_eom, new_lhs)
    return new_eom
end

"""
$(TYPEDSIGNATURES)

Transforms a higher-order differential equation system into an equivalent first-order system
by introducing additional variables. Modifies the system in place. The `time` parameter
specifies the independent variable used for differentiation.
"""
function first_order_transform!(diff_eom::DifferentialEquation, time)
    eqs′, states′ = ode_order_lowering(diff_eom.equations, time, diff_eom.harmonics)
    diff_eom.equations = eqs′
    diff_eom.harmonics = states′
    return nothing
end

"""
$(TYPEDSIGNATURES)

Helper function that performs the transformation of a higher-order differential equation system
into an equivalent first-order system. Returns a tuple of the transformed equations and the
corresponding harmonics dictionary for the new variables. Used internally by
`first_order_transform!`.
"""
function ode_order_lowering(equations, iv, harmonics)
    states = unwrap.(collect(keys(harmonics)))
    eqs = unwrap.(collect(values(equations)))

    var_order = OrderedDict{Any,Int}()
    D = Differential(iv)
    diff_eqs = empty(equations)
    diff_vars = empty(harmonics)

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
