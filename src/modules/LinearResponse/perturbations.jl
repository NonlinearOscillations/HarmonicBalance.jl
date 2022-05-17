export get_Jacobian

"""
$(SIGNATURES)

Obtain the symbolic Jacobian matrix of `eom` (either a `HarmonicEquation` or a `DifferentialEquation`).
This is the linearised left-hand side of F(u) = du/dT.

"""
function get_Jacobian(eom::HarmonicEquation)

    rearr = !HarmonicBalance.is_rearranged(eom) ? HarmonicBalance.rearrange_standard(eom) : eom
    lhs = _remove_brackets(rearr)
    vars = _remove_brackets.(eom.variables)
    
    get_Jacobian(lhs, vars)

end

" Obtain a Jacobian from a `DifferentialEquation` by first converting it into a `HarmonicEquation`. "
function get_Jacobian(diff_eom::DifferentialEquation)
    @variables T
    harmonic_eq = get_harmonic_equations(diff_eom, slow_time=T, fast_time=first(get_independent_variables(diff_eom)))
    get_Jacobian(harmonic_eq)
end


" Get the Jacobian of a set of equations `eqs` with respect to the variables `vars`. "
function get_Jacobian(eqs::Vector{Num}, vars::Vector{Num})
    length(eqs) == length(vars) || error("Jacobians are only defined for square systems!")
    M = Matrix{Num}(undef, length(vars), length(vars))

    for idx in CartesianIndices(M)
        M[idx] = expand_derivatives(d(eqs[idx[1]], vars[idx[2]]))
    end
    M
end

get_Jacobian(eqs::Vector{Equation}, vars::Vector{Num}) = get_Jacobian(Num.(getfield.(eqs, :lhs) .- getfield.(eqs, :rhs)), vars)


# the Jacobian of some equations again
function get_Jacobian_steady(eom::HarmonicEquation; differential_order=0)
    Hopf_vars = first.(getfield.(filter(x -> x.type == "Hopf", eom.variables), :symbol))
    Hopf_idx = findall(x -> any(isequal.(x, Hopf_vars)) , get_variables(eom))
    nonsingular = filter( x -> x ∉ Hopf_idx, 1:length(get_variables(eom)))

    vars_simp = Dict([var => HarmonicBalance.declare_variable(var_name(var)) for var in get_variables(eom)])
    T = get_independent_variables(eom)[1]
    J = get_Jacobian(eom.equations, d(get_variables(eom), T, differential_order))[nonsingular, nonsingular]
    
    expand_derivatives.(HarmonicBalance.substitute_all(J, vars_simp))
end

# COMPILE THIS!
function get_implicit_Jacobian(eom::HarmonicEquation)
    M = get_Jacobian_steady(eom, differential_order=0)
    Mp = get_Jacobian_steady(eom, differential_order=1)
    function J(soln::OrderedDict)
        -inv(ComplexF64.(substitute_all(Mp, soln))) * ComplexF64.(substitute_all(M, soln))
    end
    J
end


###
# STUFF BELOW IS FOR THE CORRECTED JACOBIAN METHOD
###


"""
    get_response_matrix(diff_eq::DifferentialEquation, freq::Num; order=2)

Obtain the symbolic linear response matrix of a `diff_eq` corresponding to a perturbation frequency `freq`.
This routine cannot accept a `HarmonicEquation` since there, some time-derivatives are already dropped.
`order` denotes the highest differential order to be considered.

"""
function get_response_matrix(diff_eq::DifferentialEquation, freq::Num; order=2)::Matrix
    @variables T, i
    time = get_independent_variables(diff_eq)[1]

    eom = HarmonicBalance.harmonic_ansatz(diff_eq, time)

    # replace the time-dependence of harmonic variables by slow time BUT do not drop any derivatives
    eom = HarmonicBalance.slow_flow(eom, fast_time=time, slow_time=T, degree=order+1)

    eom = HarmonicBalance.fourier_transform(eom, time)

    M = get_Jacobian(eom.equations, get_variables(eom))
    for n in 1:order
        M += (i * freq)^n * get_Jacobian(eom.equations, d(get_variables(eom), T, n))
    end
    M = substitute_all(M, [var => HarmonicBalance.declare_variable(var_name(var)) for var in get_variables(eom)])
    substitute_all(expand_derivatives.(M), i => im)
end
