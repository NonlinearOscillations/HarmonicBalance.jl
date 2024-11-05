"""
Here stability and linear response is treated with the slow-flow approximation (SFA), see Chapter 5 of JK's thesis.
Linear response always appears as a sum of Lorentzians, but is inaccurate where these are peaked far from the drive frequency.
The Jacobian is stored in the Problem object as a function that takes a solution dictionary to give the numerical Jacobian.
"""

"""
$(SIGNATURES)

Obtain the symbolic Jacobian matrix of `eom` (either a `HarmonicEquation` or a `DifferentialEquation`).
This is the linearised left-hand side of F(u) = du/dT.

"""
function get_Jacobian(eom::HarmonicEquation)
    rearr = !is_rearranged(eom) ? rearrange_standard(eom) : eom
    lhs = _remove_brackets(rearr)
    vars = _remove_brackets.(eom.variables)

    return get_Jacobian(lhs, vars)
end

" Obtain a Jacobian from a `DifferentialEquation` by first converting it into a `HarmonicEquation`. "
function get_Jacobian(diff_eom::DifferentialEquation)
    Symbolics.@variables T
    harmonic_eq = get_harmonic_equations(
        diff_eom; slow_time=T, fast_time=first(get_independent_variables(diff_eom))
    )
    return get_Jacobian(harmonic_eq)
end

" Get the Jacobian of a set of equations `eqs` with respect to the variables `vars`. "
function get_Jacobian(eqs::Vector{Num}, vars::Vector{Num})
    length(eqs) == length(vars) || error("Jacobians are only defined for square systems!")
    M = Matrix{Num}(undef, length(vars), length(vars))

    for idx in CartesianIndices(M)
        M[idx] = expand_derivatives(d(eqs[idx[1]], vars[idx[2]]))
    end
    return M
end

function get_Jacobian(eqs::Vector{Equation}, vars::Vector{Num})
    expr = Num[getfield(eq, :lhs) - getfield(eq, :rhs) for eq in eqs]
    return get_Jacobian(expr, vars)
end

"""
Code follows for an implicit treatment of the Jacobian. Usually we rearrange the linear response equations to have time-derivatives on one side.
    This may be extremely costly. Implicit evaluation means only solving the equations AFTER numerical values have been plugged in, giving
    a constant time cost per run.
"""

# for implicit evaluation, the numerical values precede the rearrangement
# for limit cycles, the zero eigenvalue causes the rearrangement to fail -> filter it out
# THIS SETS ALL DERIVATIVES TO ZERO - assumes use for steady states
function _get_J_matrix(eom::HarmonicEquation; order=0)
    order > 1 && error(
        "Cannot get a J matrix of order > 1 from the harmonic equations.\nThese are by definition missing higher derivatives",
    )

    vars_simp = Dict([var => _remove_brackets(var) for var in get_variables(eom)])
    T = get_independent_variables(eom)[1]
    J = get_Jacobian(eom.equations, d(get_variables(eom), T, order))

    return expand_derivatives.(substitute_all(J, vars_simp)) # a symbolic matrix to be compiled
end

# COMPILE THIS?
"""
$(TYPEDSIGNATURES)

Construct a function for the Jacobian of `eom` using `rules=Dict()`.

Necessary matrix inversions are only performed AFTER substituting numerical values at each call, avoiding huge symbolic operations.

Returns a function `f(soln::OrderedDict)::Matrix{ComplexF64}`.
"""
function get_implicit_Jacobian(eom::HarmonicEquation; sym_order, rules=Dict())
    J0c = compile_matrix(_get_J_matrix(eom; order=0), sym_order; rules=rules)
    J1c = compile_matrix(_get_J_matrix(eom; order=1), sym_order; rules=rules)
    m(vals::Vector) = -inv(J1c(vals)) * J0c(vals)
    return m(s::OrderedDict) = m([s[var] for var in sym_order])
end

function get_implicit_Jacobian(p::Problem, swept, fixed)
    return get_implicit_Jacobian(p.eom; sym_order=_free_symbols(p, swept), rules=fixed)
end
