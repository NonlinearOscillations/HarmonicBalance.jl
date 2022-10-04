using LinearAlgebra, Latexify
import HarmonicBalance: transform_solutions, plot
export transform_solutions, plot

""" 

    ODEProblem(
            eom::HarmonicEquation;
            fixed_parameters,
            x0::Vector,
            sweep::ParameterSweep,
            timespan::Tuple
            )

Creates an ODEProblem object used by DifferentialEquations.jl from the equations in `eom` to simulate time-evolution within `timespan`.
`fixed_parameters` must be a dictionary mapping parameters+variables to numbers (possible to use a solution index, e.g. solutions[x][y] for branch y of solution x).
If `x0` is specified, it is used as an initial condition; otherwise the values from `fixed_parameters` are used.
"""
function ODEProblem(eom::HarmonicEquation, fixed_parameters; sweep::ParameterSweep=ParameterSweep(), x0::Vector=[], timespan::Tuple, perturb_initial=0.0, kwargs...)

    if !is_rearranged(eom) # check if time-derivatives of the variable are on the right hand side
        eom = HarmonicBalance.rearrange_standard(eom)
    end

    fixed = Dict(fixed_parameters)
    # substitute fixed parameters
    fixed = HarmonicBalance.filter_duplicate_parameters(sweep, fixed)
    p_values = [fixed[p] for p in keys(fixed)]
    subeqs = substitute_all(Num.(getfield.(eom.equations, :lhs)), Dict(zip(keys(fixed), p_values)))
    vars = get_variables(eom)

    # substitute the harmonic variables
    eqs(v) = [substitute(eq, Dict(zip(vars, v))) for eq in subeqs]
     # substitute  sweep parameters
    eqs(v,T) = [substitute(eq,Dict(zip(keys(sweep), [sweep[p](T) for p in keys(sweep)]))) for eq in eqs(v)]
    
    function f!(du,u,p,T)
        state = eqs(u,T)
        for j in 1:length(vars)
            du[j] = state[j]
        end
    end

    # the initial condition is x0 if specified, taken from fixed_parameters otherwise
    initial = isempty(x0) ? real.(collect(values(fixed_parameters))[1:length(vars)]) * (1-perturb_initial) : x0

    return DifferentialEquations.ODEProblem(f!, initial, timespan; kwargs...)
end


"""
$(TYPEDSIGNATURES)

Numerically investigate the stability of a solution `soln` of `eom` within `timespan`.
The initial condition is displaced by `perturb_initial`.

Return `true` the solution evolves within `tol` of the initial value (interpreted as stable).

"""
function is_stable(soln::StateDict, eom::HarmonicEquation; timespan, tol=1E-1, perturb_initial=1E-3)
    problem = ODEProblem(eom, steady_solution=soln, timespan=timespan)
    solution = solve(problem)
    dist = norm(solution[end] - solution[1]) / (norm(solution[end]) + norm(solution[1]))
    !is_real(solution[end]) || !is_real(solution[1]) ? error("the solution is complex!") : dist < tol
end


transform_solutions(soln::OrdinaryDiffEq.ODECompositeSolution, f::String, harm_eq::HarmonicEquation) = transform_solutions(soln.u, f, harm_eq)
transform_solutions(s::OrdinaryDiffEq.ODECompositeSolution, funcs::Vector{String}, he::HarmonicEquation) = [transform_solutions(s, f, he) for f in funcs]



"""
    plot(soln::ODECompositeSolution, f::String, harm_eq::HarmonicEquation; kwargs...)

Plot a function `f` of a time-dependent solution `soln` of `harm_eq`.

## As a function of time

    plot(soln::ODECompositeSolution, f::String, harm_eq::HarmonicEquation; kwargs...)

`f` is parsed by Symbolics.jl

## parametric plots
    plot(soln::ODECompositeSolution, f::Vector{String}, harm_eq::HarmonicEquation; kwargs...)

Parametric plot of f[1] against f[2]

"""
function plot(soln::OrdinaryDiffEq.ODECompositeSolution, funcs, harm_eq::HarmonicEquation; kwargs...)
    if funcs isa String || length(funcs) == 1
        plot(soln.t, transform_solutions(soln, funcs, harm_eq); legend=false, xlabel="time", ylabel=latexify(funcs), HarmonicBalance._set_Plots_default..., kwargs...)
    elseif length(funcs) == 2 # plot of func vs func
        plot(transform_solutions(soln, funcs, harm_eq)...; legend=false, xlabel=latexify(funcs[1]), ylabel=latexify(funcs[2]), HarmonicBalance._set_Plots_default..., kwargs...)
    else
        error("Invalid plotting argument: ", funcs)
    end
end