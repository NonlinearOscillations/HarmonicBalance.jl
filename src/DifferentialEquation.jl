"""
$(TYPEDEF)

Holds differential equation(s) of motion and a set of harmonics to expand each variable.
This is the primary input for `HarmonicBalance.jl` ; after inputting the equations, the harmonics
    ansatz needs to be specified using `add_harmonic!`.

# Fields
$(TYPEDFIELDS)

## Example
```julia-repl
julia> @variables t, x(t), y(t), ω0, ω, F, k;

# equivalent ways to enter the simple harmonic oscillator
julia> DifferentialEquation(d(x,t,2) + ω0^2 * x - F * cos(ω*t), x);
julia> DifferentialEquation(d(x,t,2) + ω0^2 * x ~ F * cos(ω*t), x);

# two coupled oscillators, one of them driven
julia> DifferentialEquation([d(x,t,2) + ω0^2 * x - k*y, d(y,t,2) + ω0^2 * y - k*x] .~ [F * cos(ω*t), 0], [x,y]);
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
                    "The equation is of the form $(typerhs)~$(typelhs) is not supported. Commenly one forgot to broadcast the equation symbol `~`.",
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
                "The variables are of type $(typeof(vars)). Whereas you gave one equation of type $(typerhs)~$(typelhs). Commenly one forgot to broadcast the equation symbol `~`.",
            ),
        )
    end
    DifferentialEquation(lhs::Num, var::Num) = DifferentialEquation([lhs ~ Int(0)], [var])
end

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

"""
$(TYPEDSIGNATURES)
Add the harmonic `ω` to the harmonic ansatz used to expand the variable `var` in `diff_eom`.

## Example

# define the simple harmonic oscillator and specify that x(t) oscillates with frequency ω
```julia-repl
julia> @variables t, x(t), y(t), ω0, ω, F, k;
julia> diff_eq = DifferentialEquation(d(x,t,2) + ω0^2 * x ~ F * cos(ω*t), x);
julia> add_harmonic!(diff_eq, x, ω) # expand x using ω

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

ExprUtils.is_harmonic(diff_eom::DifferentialEquation, t::Num)::Bool =
    all([is_harmonic(eq, t) for eq in values(diff_eom.equations)])

"Pretty printing of the newly defined types"
function show_fields(object)
    for field in fieldnames(typeof(object)) # display every field
        display(string(field))
        display(getfield(object, field))
    end
end

"""
$(TYPEDSIGNATURES)
Return the independent dependent variables of `diff_eom`.
"""
function get_independent_variables(diff_eom::DifferentialEquation)
    return Num.(flatten(unique([x.val.arguments for x in keys(diff_eom.equations)])))
end

Base.show(eom::DifferentialEquation) = show_fields(eom)
