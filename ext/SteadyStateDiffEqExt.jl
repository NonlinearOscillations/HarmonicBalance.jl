module SteadyStateDiffEqExt

export steady_state_sweep
import HarmonicBalance: steady_state_sweep
using SteadyStateDiffEq: solve, NonlinearProblem, SteadyStateProblem, DynamicSS, remake
using LinearAlgebra: norm, eigvals
using SteadyStateDiffEq.SciMLBase.SciMLStructures: isscimlstructure

function steady_state_sweep(
        prob::SteadyStateProblem, alg::DynamicSS;
        varied::Pair, kwargs...)
    varied_idx, sweep_range = varied

    # if p is dual number (AD) result is dual number
    result = [similar(prob.u0) for _ in sweep_range]

    foreach(pairs(sweep_range)) do (i, value)
        u0 = i == 1 ? [0.0, 0.0] : result[i - 1]
        # make type-stable: FD.Dual or Float64
        parameters = get_new_parameter(prob, varied_idx, value)
        sol = solve(remake(prob, p = parameters, u0 = u0), alg; kwargs...)
        result[i] = sol.u
    end
    return result
end

function steady_state_sweep(
        prob_np::NonlinearProblem, prob_ss::SteadyStateProblem,
        alg_np, alg_ss::DynamicSS;
        varied::Pair, kwargs...)
    varied_idx, sweep_range = varied
    # if p is dual number (AD) result is dual number
    result = [similar(prob_np.u0) for _ in sweep_range]

    foreach(pairs(sweep_range)) do (i, value)
        u0 = i == 1 ? Base.zeros(length(prob_np.u0)) : result[i - 1]
        # make type-stable: FD.Dual or Float64
        parameters = get_new_parameter(prob_np, varied_idx, value)
        sol_nn = solve(remake(prob_np, p = parameters, u0 = u0), alg_np; kwargs...)

        # last argument is time but does not matter
        zeros = norm(prob_np.f.f.f.f.f_oop(sol_nn.u, parameters, 0))
        jac = prob_np.f.jac.f.f.f_oop(sol_nn.u, parameters, 0)
        eigval = jac isa Vector ? jac : eigvals(jac) # eigvals favourable supports FD.Dual

        if !isapprox(zeros, 0, atol = 1e-5) || any(λ -> λ > 0, real.(eigval))
            sol_ss = solve(remake(prob_ss, p = parameters, u0 = u0),
                alg_ss, abstol = 1e-5, reltol = 1e-5)
            result[i] = sol_ss.u
        else
            result[i] = sol_nn.u
        end
    end

    return result
end

function get_new_parameter(prob, varied_idx, value)
    # make type-stable: FD.Dual or Float64
    if hasfield(typeof(prob.p), :tunable)
        rest = [prob.p.discrete, prob.p.nonnumeric, prob.p.dependent, prob.p.constant]

        all(isempty.(rest)) || error("Only tunable parameters are supported")
        length(prob.p.tunable) == 1 ||
            error("The type of the parameters should be uniform")

        old_parameters_values = prob.p.tunable[1]
        parameter_values = eltype(old_parameters_values)[i == varied_idx ? value : x
                                          for (i, x) in enumerate(old_parameters_values)]
        parameters = replace(Tunable(), prob.p, parameter_values)
    else
        parameters = eltype(prob.p)[i == varied_idx ? value : x
                                    for (i, x) in enumerate(prob.p)]
    end
    return parameters
end

end # module
