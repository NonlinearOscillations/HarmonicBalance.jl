using OrdinaryDiffEqTsit5: ODEProblem, solve, ODESolution

"""
    $(TYPEDSIGNATURES)

Creates an ODEProblem object used by OrdinaryDiffEqTsit5.jl from the equations in `eom` to simulate time-evolution within `timespan`.
`fixed_parameters` must be a dictionary mapping parameters+variables to numbers (possible to use a solution index, e.g. solutions[x][y] for branch y of solution x).
If `u0` is specified, it is used as an initial condition; otherwise the values from `fixed_parameters` are used.
"""
function OrdinaryDiffEqTsit5.ODEProblem(
    eom::HarmonicEquation,
    fixed_parameters;
    sweep::AdiabaticSweep=AdiabaticSweep(),
    u0::Vector=[],
    timespan::Tuple,
    perturb_initial=0.0,
    kwargs...,
)
    if !is_rearranged(eom) # check if time-derivatives of the variable are on the right hand side
        eom = rearrange_standard(eom)
    end

    fixed = Dict(fixed_parameters)
    # substitute fixed parameters
    fixed = filter_duplicate_parameters(sweep, fixed)
    p_values = [fixed[p] for p in keys(fixed)]
    subeqs = substitute_all(
        Num.(getfield.(eom.equations, :lhs)), Dict(zip(keys(fixed), p_values))
    )
    vars = get_variables(eom)

    # substitute the harmonic variables
    eqs(v) = [substitute(eq, Dict(zip(vars, v))) for eq in subeqs]
    # substitute  sweep parameters
    function eqs(v, T)
        return real.(
            unwrap.([
                substitute(eq, Dict(zip(keys(sweep), [sweep[p](T) for p in keys(sweep)])))
                for eq in eqs(v)
            ])
        )
    end

    function f!(du, u, p, T) # in-place
        state = eqs(u, T)
        for j in 1:length(vars)
            du[j] = state[j]
        end
    end

    # the initial condition is u0 if specified, taken from fixed_parameters otherwise
    initial = if isempty(u0)
        real.(collect(values(fixed_parameters))[1:length(vars)]) * (1 - perturb_initial)
    else
        u0
    end

    return ODEProblem(f!, initial, timespan; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Numerically investigate the stability of a solution `soln` of `eom` within `timespan`.
The initial condition is displaced by `perturb_initial`.

Return `true` the solution evolves within `tol` of the initial value (interpreted as stable).

"""
function HarmonicBalance.is_stable(
    steady_solution::StateDict,
    eom::HarmonicEquation;
    timespan,
    tol=1e-1,
    perturb_initial=1e-3,
)
    problem = ODEProblem(eom; steady_solution, timespan)
    solution = solve(problem)
    dist = norm(solution[end] - solution[1]) / (norm(solution[end]) + norm(solution[1]))
    return if !is_real(solution[end]) || !is_real(solution[1])
        error("the solution is complex!")
    else
        dist < tol
    end
end

function HarmonicBalance.transform_solutions(
    soln::ODESolution, f::String, harm_eq::HarmonicEquation
)
    return transform_solutions(soln.u, f, harm_eq)
end
function HarmonicBalance.transform_solutions(
    s::ODESolution, funcs::Vector{String}, harm_eq::HarmonicEquation
)
    return [transform_solutions(s, f, harm_eq) for f in funcs]
end
