module SteadyStateDiffEqExt

export steady_state_sweep
using HarmonicBalance: HarmonicBalance, steady_state_sweep

using SteadyStateDiffEq: solve, NonlinearProblem, SteadyStateProblem, DynamicSS, remake
using LinearAlgebra: norm, eigvals
using SteadyStateDiffEq.SciMLBase.SciMLStructures: Tunable, replace
# using NonlinearSolve.NonlinearSolveBase.Utils: evaluate_f
# evaluate_f(prob_np,[0])

"""
    steady_state_sweep(prob::SteadyStateProblem, alg::DynamicSS; varied::Pair, kwargs...)

Sweeps through a range of parameter values using a dynamic steady state solver `DynamicSS`
of the `SteadyStateDiffEq.jl` package. Given a steady state problem and a parameter to vary,
computes the steady state solution for each value in the sweep range. The solutions are
returned as a vector where each element corresponds to the steady state found at that
parameter value.
"""
function HarmonicBalance.steady_state_sweep(
    prob::SteadyStateProblem, alg::DynamicSS; varied::Pair, kwargs...
)
    varied_idx, sweep_range = varied

    # if p is dual number (AD) result is dual number
    result = [similar(prob.u0) for _ in sweep_range]

    foreach(pairs(sweep_range)) do (i, value)
        u0 = i == 1 ? [0.0, 0.0] : result[i - 1]
        # make type-stable: FD.Dual or Float
        p = get_new_parameters(prob, varied_idx, value)
        sol = solve(remake(prob; p, u0), alg; kwargs...)
        result[i] = sol.u
    end
    return result
end

"""
    steady_state_sweep(prob_np::NonlinearProblem, prob_ss::SteadyStateProblem,
                      alg_np, alg_ss::DynamicSS; varied::Pair, kwargs...)

Performs a parameter sweep by combining nonlinear root `alg_np` and steady state solvers `alg_ss`.
For each parameter value, it first attempts a direct nonlinear root solver and checks its
stability. If the solution is unstable or not found, it switches to a dynamic steady state solver.
This hybrid approach is much faster then only using a steady state solver.
"""
function HarmonicBalance.steady_state_sweep(
    prob_np::NonlinearProblem,
    prob_ss::SteadyStateProblem,
    alg_np,
    alg_ss::DynamicSS;
    varied::Pair,
    kwargs...,
)
    varied_idx, sweep_range = varied
    # if p is dual number (AD) result is dual number
    result = [similar(prob_np.u0) for _ in sweep_range]

    foreach(pairs(sweep_range)) do (i, value)
        u0 = i == 1 ? Base.zeros(length(prob_np.u0)) : result[i - 1]
        # make type-stable: FD.Dual or Float
        p = get_new_parameters(prob_np, varied_idx, value)
        sol_nn = solve(remake(prob_np; p, u0), alg_np; kwargs...)

        param_val = tunable_parameters(p)
        zeros = norm(prob_ss.f(sol_nn.u, param_val, NaN)) # todo check NaN
        jac = prob_ss.f.jac(sol_nn.u, param_val, NaN)
        eigval = jac isa AbstractVector ? jac : eigvals(jac) # eigvals favourable supports FD.Dual

        if !isapprox(zeros, 0; atol=1e-5) || any(λ -> λ > 0, real.(eigval))
            sol_ss = solve(remake(prob_ss; p, u0), alg_ss; abstol=1e-5, reltol=1e-5)
            result[i] = sol_ss.u
        else
            result[i] = sol_nn.u
        end
    end

    return result
end

function tunable_parameters(param)
    return hasfield(typeof(param), :tunable) ? param.tunable : param
end

function get_new_parameters(prob, varied_idx, value)
    # make type-stable: FD.Dual or Float
    if hasfield(typeof(prob.p), :tunable)
        rest = map(filter(x -> x != :tunable, propertynames(prob.p))) do prop
            getproperty(prob.p, prop)
        end

        all(isempty.(rest)) || error("Only tunable parameters are supported")

        old_parameters_values = prob.p.tunable
        parameter_values = eltype(old_parameters_values)[
            i == varied_idx ? value : x for (i, x) in enumerate(old_parameters_values)
        ]
        parameters = replace(Tunable(), prob.p, parameter_values)
    else
        parameters = eltype(prob.p)[
            i == varied_idx ? value : x for (i, x) in enumerate(prob.p)
        ]
    end
    return parameters
end

end # module
