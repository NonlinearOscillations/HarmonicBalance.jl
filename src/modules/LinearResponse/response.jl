"""
    get_response_matrix(diff_eq::DifferentialEquation, freq::Num; order=2)

Obtain the symbolic linear response matrix of a `diff_eq` corresponding to a perturbation frequency `freq`.
This routine cannot accept a `HarmonicEquation` since there, some time-derivatives are already dropped.
`order` denotes the highest differential order to be considered.

"""
function get_response_matrix(diff_eq::DifferentialEquation, freq::Num; order=2)::Matrix
    @variables T, i
    time = get_independent_variables(diff_eq)[1]

    eom = harmonic_ansatz(diff_eq, time)

    # replace the time-dependence of harmonic variables by slow time BUT do not drop any derivatives
    eom = slow_flow(eom; fast_time=time, slow_time=T, degree=order + 1)

    eom = fourier_transform(eom, time)

    # get the response matrix by summing the orders
    M = get_Jacobian(eom.equations, get_variables(eom))
    for n in 1:order
        M += (i * freq)^n * get_Jacobian(eom.equations, d(get_variables(eom), T, n))
    end
    M = substitute_all(
        M, [var => declare_variable(var_name(var)) for var in get_variables(eom)]
    )
    return substitute_all(expand_derivatives.(M), i => im)
end

"Get the response matrix corresponding to `res`.
Any substitution rules not specified in `res` can be supplied in `rules`."
function ResponseMatrix(res::Result; rules=Dict())

    # get the symbolic response matrix
    @variables Δ
    M = get_response_matrix(res.problem.eom.natural_equation, Δ; order=2)
    M = substitute_all(M, merge(res.fixed_parameters, rules))
    symbols = _free_symbols(res)
    compiled_M = [build_function(el, cat(symbols, [Δ]; dims=1)) for el in M]

    return ResponseMatrix(eval.(compiled_M), symbols, res.problem.eom.variables)
end

"""Evaluate the response matrix `resp` for the steady state `s` at (lab-frame) frequency `Ω`."""
function evaluate_response_matrix(resp::ResponseMatrix, s::StateDict, Ω)
    values = cat([s[var] for var in resp.symbols], [Ω]; dims=1)
    f = resp.matrix
    return [Base.invokelatest(el, values) for el in f]
end

## THIS NEEDS REVISING
function _evaluate_response_vector(rmat::ResponseMatrix, s::StateDict, Ω)
    m = evaluate_response_matrix(rmat, s, Ω)
    force_pert = cat([[1.0, 1.0 * im] for n in 1:(size(m)[1] / 2)]...; dims=1)
    return inv(m) * force_pert
end

"""
$(TYPEDSIGNATURES)

For `rmat` and a solution dictionary `s`,
calculate the total response to a perturbative force at frequency `Ω`.

"""
function get_response(rmat::ResponseMatrix, s::StateDict, Ω)
    resp = 0

    # uv-type
    for pair in _get_uv_pairs(rmat.variables)
        u, v = rmat.variables[pair]
        this_ω = unwrap(substitute_all(u.ω, s))
        uv1 = _evaluate_response_vector(rmat, s, Ω - this_ω)[pair]
        uv2 = _evaluate_response_vector(rmat, s, -Ω + this_ω)[pair]
        resp += sqrt(_plusamp(uv1)^2 + _minusamp(uv2)^2)
    end

    # a-type variables
    for a_idx in _get_as(rmat.variables)
        a = rmat.variables[a_idx]
        uv1 = _evaluate_response_vector(rmat, s, Ω)[a]
        uv2 = _evaluate_response_vector(rmat, s, -Ω)[a]
        resp += sqrt(_plusamp(uv1)^2 + _minusamp(uv2)^2)
    end
    return resp
end

# formulas to obtain up- and down- converted frequency components when going from the
# rotating frame into the lab frame
_plusamp(uv) = norm(uv)^2 - 2 * (imag(uv[1]) * real(uv[2]) - real(uv[1]) * imag(uv[2]))
_minusamp(uv) = norm(uv)^2 + 2 * (imag(uv[1]) * real(uv[2]) - real(uv[1]) * imag(uv[2]))
