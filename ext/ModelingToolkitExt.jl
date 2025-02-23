module ModelingToolkitExt

using DocStringExtensions

export ODESystem, ODEProblem, SteadyStateProblem, NonlinearProblem

using HarmonicBalance:
    HarmonicEquation,
    is_rearranged,
    rearrange_standard,
    get_variables,
    DifferentialEquation,
    get_independent_variables
using HarmonicBalance.KrylovBogoliubov:
    rearrange_standard!, is_rearranged_standard, first_order_transform!
using Symbolics:
    simplify, Equation, substitute, Num, @variables, expand, unwrap, arguments, wrap
using ModelingToolkit:
    ModelingToolkit,
    ODESystem,
    ODEProblem,
    NonlinearProblem,
    SteadyStateProblem,
    varmap_to_vars,
    parameters,
    @parameters,
    @mtkbuild,
    @independent_variables

swapsides(eq::Equation) = Equation(eq.rhs, eq.lhs)

function declare_parameter(var::Num)
    var_sym = Symbol(var)
    new_var = @parameters $var_sym
    @eval($(var_sym) = first($new_var)) # store the variable under "name" in this namespace
    return eval(var_sym)
end

@doc """
$(TYPEDSIGNATURES)

Creates and ModelingToolkit.ODESystem from a HarmonicEquation.

### Example
```julia
using ModelingToolkit

@variables α ω ω0 F γ t x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω) #
harmonic_eq = get_harmonic_equations(diff_eq)

sys = ODESystem(harmonic_eq)
param = (α => 1.0, ω0 => 1.1, F => 0.01, γ => 0.01, ω => 1.1)
ODEProblem(sys, [1.0, 0.0], (0, 100), param)
```
"""
function ModelingToolkit.ODESystem(eom::HarmonicEquation)
    if !is_rearranged(eom) # check if time-derivatives of the variable are on the right hand side
        eom = rearrange_standard(eom)
    end

    vars = get_variables(eom)
    slow_times = arguments.(unwrap.(vars))
    @assert all(isone.(length.(slow_times))) "Only one argument for the variables are allowed."
    slow_time = unique(first.(slow_times))
    @assert isone(length(slow_time)) "The argument of the variables are not the same."
    slow_time_ivp = @eval @independent_variables $(Symbol(first(slow_time)))

    # we have the replace the param made with @variables with the ones made with @parameters
    par_names = declare_parameter.(eom.parameters)

    eqs = deepcopy(eom.equations)
    eqs = swapsides.(eqs)
    eqs = simplify.(expand.(eqs))
    eqs = substitute(eqs, Dict(zip(eom.parameters, par_names)))

    # ∨ mtk v9 need @mtkbuild
    @mtkbuild sys = ODESystem(eqs, first(slow_time_ivp), vars, par_names)
    return sys
end

@doc """
$(TYPEDSIGNATURES)

Creates and ModelingToolkit.ODESystem from a DifferentialEquation.

### Example
```julia
using ModelingToolkit

@variables α ω ω0 F γ t x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
sys = ODESystem(diff_eq)

param = (α => 1.0, ω0 => 1.1, F => 0.01, γ => 0.01, ω => 1.1)

ODEProblem(sys, [1.0, 0.0], (0, 100), param)
```
"""
function ModelingToolkit.ODESystem(diff_eq::DifferentialEquation)
    diff_eq = deepcopy(diff_eq)
    if !is_rearranged_standard(diff_eq)
        rearrange_standard!(diff_eq)
    end

    times = get_independent_variables(diff_eq)
    @assert isone(length(times)) "Only one independent variable allowed."
    iv = first(@eval @independent_variables $(Symbol(first(times))))

    first_order_transform!(diff_eq, iv)

    eqs = collect(values(diff_eq.equations))
    vars = get_variables(diff_eq)

    diff_eq_sym = collect(Iterators.flatten(get_variables.(eqs)))
    param_undeclared = setdiff(setdiff(wrap.(diff_eq_sym), vars), iv)
    params = declare_parameter.(param_undeclared)

    eqs = substitute(eqs, Dict(zip(param_undeclared, params)))

    @mtkbuild sys = ODESystem(eqs, first(iv), vars, params)

    return sys
end

@doc """
$(TYPEDSIGNATURES)

Creates and ModelingToolkit.ODEProblem from a DifferentialEquation or HarmonicEquation.

### Example
```julia
using ModelingToolkit, StaticArrays

@variables α ω ω0 F γ t x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω) #
harmonic_eq = get_harmonic_equations(diff_eq)

# in place (most performant for large systems)
ODEProblem(harmonic_eq, [1.0, 0.0], (0, 100), param)

# out of place (most performant for small systems with StaticArrays)
ODEProblem(
    harmonic_eq, [1.0, 0.0], (0, 100), param;
    in_place=false, u0_constructor=x -> SVector(x...)
)
```

"""
function ModelingToolkit.ODEProblem(
    eom::Union{HarmonicEquation,DifferentialEquation},
    u0,
    tspan::Tuple,
    p::AbstractDict;
    in_place=true,
    kwargs...,
)
    sys = ODESystem(eom)
    param = varmap_to_vars(p, parameters(sys))
    if !in_place # out-of-place
        prob = ODEProblem{false}(sys, u0, tspan, param; jac=true, kwargs...)
    else # in-place
        prob = ODEProblem{true}(sys, u0, tspan, param; jac=true, kwargs...)
    end # compute jacobian for performance
    return prob
end

@doc """
$(TYPEDSIGNATURES)

Creates and ModelingToolkit.NonlinearProblem from a HarmonicEquation.

### Example
```julia
using ModelingToolkit, StaticArrays

@variables α ω ω0 F γ t x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω) #
harmonic_eq = get_harmonic_equations(diff_eq)


NonlinearProblem(harmonic_eq, [1.0, 0.0], param)
```
"""
function ModelingToolkit.NonlinearProblem(
    eom::HarmonicEquation, u0, p::AbstractDict; in_place=true, kwargs...
)
    ss_prob = SteadyStateProblem(eom, u0, p::AbstractDict; in_place, kwargs...)
    return NonlinearProblem(ss_prob)
end

@doc """
$(TYPEDSIGNATURES)

Creates and ModelingToolkit.SteadyStateProblem from a HarmonicEquation.

### Example
```julia
using ModelingToolkit, StaticArrays

@variables α ω ω0 F γ t x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω) #
harmonic_eq = get_harmonic_equations(diff_eq)


SteadyStateProblem(harmonic_eq, [1.0, 0.0], param)
```
"""
function ModelingToolkit.SteadyStateProblem(
    eom::HarmonicEquation, u0, p::AbstractDict; in_place=true, kwargs...
)
    sys = ODESystem(eom)
    param = varmap_to_vars(p, parameters(sys))
    if !in_place # out-of-place
        prob = SteadyStateProblem{false}(sys, u0, param; jac=true, kwargs...)
    else # in-place
        prob = SteadyStateProblem{true}(sys, u0, param; jac=true, kwargs...)
    end # compute jacobian for performance
    return prob
end

end # module
