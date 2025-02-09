"""
Compile the Jacobian from `prob`, inserting `fixed_parameters`.
Returns a function that takes a vector of variables and `swept_parameters` to give the Jacobian.
The order of the vector is first the variables, then the swept parameters.
For LC variables called "hopf", the Jacobian is already compiled.
For Type stability the function is wrapped in a FunctionWrapper.
"""
function _compile_Jacobian(
    eom::HarmonicEquation, soltype::DataType, swept::OrderedDict, fixed::OrderedDict
)::JacobianFunction(soltype)
    if "Hopf" ∈ getfield.(eom.variables, :type)
        compiled_J = eom.jacobian
    elseif !hasnan(eom.jacobian)
        compiled_J = compile_matrix(eom.jacobian, _free_symbols(eom, swept); rules=fixed)
    else
        compiled_J = get_implicit_Jacobian(
            eom; sym_order=_free_symbols(eom, swept), rules=fixed
        )
    end
    return JacobianFunction(soltype)(compiled_J)
end

"""
Take a matrix containing symbolic variables `variables` and keys of `fixed_parameters`.
Substitute the values according to `fixed_parameters` and compile into a function that takes
numerical arguments in the order set in `variables`.
"""
function compile_matrix(
    mat::Matrix{Num}, variables::Vector{Num}; rules=Dict()
)::Symbolics.RuntimeGeneratedFunction
    J = substitute_all.(mat, Ref(rules)) # Ref makes sure only mat is broadcasted
    jacfunc = Symbolics.build_function(J, variables; expression=Val(false))
    return jacfunc isa Tuple ? jacfunc[1] : jacfunc
end

"""
The Jacobian is stored in the Problem object as a
function that takes a solution dictionary to give the numerical Jacobian.
"""

"""
$(SIGNATURES)

Obtain the symbolic Jacobian matrix of `eom`.
This is the linearised left-hand side of F(u) = du/dT.
"""
function get_Jacobian(eom::HarmonicEquation)::Matrix{Num}
    rearr = !is_rearranged(eom) ? rearrange_standard(eom) : eom
    lhs = _remove_brackets(rearr)
    vars = _remove_brackets.(eom.variables)

    return get_Jacobian(lhs, vars)
end

function add_jacobian!(eom::HarmonicEquation)
    return eom.jacobian .= get_Jacobian(eom)
end

" Obtain a Jacobian from a `DifferentialEquation` by first converting it into a `HarmonicEquation`. "
function get_Jacobian(diff_eom::DifferentialEquation)::Matrix{Num}
    Symbolics.@variables T
    harmonic_eq = get_harmonic_equations(
        diff_eom; slow_time=T, fast_time=first(get_independent_variables(diff_eom))
    )
    return get_Jacobian(harmonic_eq)
end

" Get the Jacobian of a set of equations `eqs` with respect to the variables `vars`. "
function get_Jacobian(eqs::Vector{Num}, vars::Vector{Num})::Matrix{Num}
    length(eqs) == length(vars) || error("Jacobians are only defined for square systems!")
    M = Matrix{Num}(undef, length(vars), length(vars))

    for idx in CartesianIndices(M)
        M[idx] = expand_derivatives(d(eqs[idx[1]], vars[idx[2]]))
    end
    return M
end

function get_Jacobian(eqs::Vector{Equation}, vars::Vector{Num})::Matrix{Num}
    expr = Num[getfield(eq, :lhs) - getfield(eq, :rhs) for eq in eqs]
    return get_Jacobian(expr, vars)
end

"""
Code follows for an implicit treatment of the Jacobian. Usually we rearrange the linear response
equations to have time-derivatives on one side. This may be extremely costly. Implicit evaluation
means only solving the equations AFTER numerical values have been plugged in, giving a constant
time cost per run.
"""

# for implicit evaluation, the numerical values precede the rearrangement
# for limit cycles, the zero eigenvalue causes the rearrangement to fail -> filter it out
# THIS SETS ALL DERIVATIVES TO ZERO - assumes use for steady states
function _get_J_matrix(eom::HarmonicEquation; order=0)
    order > 1 && error("Cannot get a J matrix of order > 1 from the harmonic equations.\n
                       These are by definition missing higher derivatives")

    vars_simp = Dict([var => _remove_brackets(var) for var in get_variables(eom)])
    T = get_independent_variables(eom)[1]
    J = get_Jacobian(eom.equations, d(get_variables(eom), T, order))

    return expand_derivatives.(substitute_all(J, vars_simp)) # a symbolic matrix to be compiled
end

# TODO COMPILE THIS?
"""
$(TYPEDSIGNATURES)

Construct a function for the Jacobian of `eom` using `rules=Dict()`.

Necessary matrix inversions are only performed AFTER substituting numerical values at each call,
avoiding huge symbolic operations.

Returns a function `f(soln::OrderedDict{Num,T})::Matrix{T}`.
"""
function get_implicit_Jacobian(eom::HarmonicEquation; sym_order, rules=Dict())
    J0c = compile_matrix(_get_J_matrix(eom; order=0), sym_order; rules)
    J1c = compile_matrix(_get_J_matrix(eom; order=1), sym_order; rules)
    jacfunc(vals::Vector) = -inv(real.(J1c(vals))) * J0c(vals)
    return jacfunc
end

function get_implicit_Jacobian(p::Problem)
    return get_implicit_Jacobian(
        p.eom; sym_order=_free_symbols(p), rules=p.fixed_parameters
    )
end

function dummy_symbolic_Jacobian(n::Int)::Matrix{Num}
    return Num.(float.(collect(LinearAlgebra.I(n))) .* NaN)
end
