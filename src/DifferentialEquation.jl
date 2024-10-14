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
