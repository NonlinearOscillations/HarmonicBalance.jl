using LinearAlgebra


""" 

    ODEProblem(
            eom::HarmonicEquation;
            fixed_parameters::ParameterList,
            x0::Vector,
            steady_solution::Dict
            sweep::ParameterSweep,
            timespan::Tuple
            )

Creates an ODEProblem object used by DifferentialEquations.jl from the equations in `eom` to simulate time-evolution within `timespan`.
To manually input parameter values and initial conditions, use the keywords `fixed_parameters` and `x0`.
To start the evolution from a steady-state solution, use `steady_solution`.

"""
function ODEProblem(eom::HarmonicEquation, fixed_parameters; sweep::ParameterSweep=ParameterSweep(), x0::Vector, timespan::Tuple)

    if !is_rearranged(eom) # check if time-derivatives of the variable are on the right hand side
        eom = HarmonicBalance.rearrange_standard(eom)
    end

    # substitute fixed parameters
    fixed_parameters = HarmonicBalance.filter_duplicate_parameters(sweep, fixed_parameters)
    p_values = [fixed_parameters[p] for p in keys(fixed_parameters)]
    subeqs = substitute_all(Num.(getfield.(eom.equations, :lhs)), Dict(zip(keys(fixed_parameters), p_values)))
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

    return DifferentialEquations.ODEProblem(f!,x0,timespan)
end

# evolving from a steady-state solution found with homotopy continuation
function ODEProblem(eom::HarmonicEquation; steady_solution::StateDict, sweep=ParameterSweep(), timespan, perturb_initial=0)
    vars = [HarmonicBalance.declare_variable(v) for v in (HarmonicBalance.var_name.(get_variables(eom)))]
    initial = real.([steady_solution[v] for v in vars]) * (1-perturb_initial)
    ODEProblem(eom, steady_solution, sweep=sweep, x0=initial, timespan=timespan)
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