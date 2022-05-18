export add_Hopf!
export Hopf_variables
export replace_Hopf_variable!

using HarmonicBalance: is_rearranged, rearrange_standard
using HarmonicBalance: LinearResponse.get_implicit_Jacobian


"Add a Hopf bifurcation to the system."
function add_Hopf!(eom::DifferentialEquation, Δω::Num)
    for var in get_variables(eom), ω in eom.harmonics[var]
        add_harmonic!(eom, var, ω + Δω)
        add_harmonic!(eom, var, ω - Δω)
    end
end


"Successively adds n Hopf bifurcations at the same frequency Δω"
add_Hopf!(eom::DifferentialEquation; frequency::Num, multiplicity::Int) = [add_Hopf!(eom, frequency) for k in 1:multiplicity]


"""
    Hopf_variables(eom, freq)
"""
Hopf_variables(eom::HarmonicEquation, freq::Num) = filter(x -> any(isequal.(freq, get_all_terms(x.ω))), eom.variables)


"""
    Obtain the Jacobian of `eom` with a free (Hopf) variable `free_var`.
    `free_var` marks the variable which is free due to U(1) symmetry. Its entry in the Jacobian matrix must be discarded
    (the discarded eigenvalue is 0, corresponding to a free phase.)
"""
function _Hopf_Jacobian(eom::HarmonicEquation, free_var::HarmonicVariable)
    eom_Jac = rearrange_standard(eom)
    free_idx = findall(x -> isequal(x, free_var), eom_Jac.variables)  
    deleteat!(eom_Jac.equations, free_idx)
    deleteat!(eom_Jac.variables, free_idx)

    # the free variable can be fixed to zero, this is also done in the corresponding HarmonicEquation later
    eom_Jac = substitute_all(eom_Jac, free_var => 0)
    J = get_implicit_Jacobian(eom_Jac)
end


"""
    Construct a `Problem` in the case where U(1) symmetry is present
    due to having added a Hopf variable with frequency `Hopf_ω`.
"""
function _Hopf_Problem(eom::HarmonicEquation, Hopf_ω::Num)

    isempty(Hopf_variables(eom, Hopf_ω)) ? error("No Hopf variables found!") : nothing
    !any(isequal.(eom.parameters, Hopf_ω)) ? error(Hopf_ω, " is not a parameter of the harmonic equation!") : nothing 

    free_var = Hopf_variables(eom, Hopf_ω)[end] # eliminate one of the Cartesian variables, it does not matter which
    
    # get the Hopf Jacobian before altering anything - this is the usual Jacobian but the entries corresponding
    # to the free variable are removed
    J = _Hopf_Jacobian(eom, free_var)
    
    new_symbol = HarmonicBalance.declare_variable(string(Hopf_ω), first(get_independent_variables(eom)))
    rules = Dict(Hopf_ω => new_symbol, free_var.symbol => Num(0))
    eom.equations = expand_derivatives.(substitute_all(eom.equations, rules))
    eom.parameters = setdiff(eom.parameters, [Hopf_ω]) # Hopf_ω is now NOT a parameter anymore
    
    free_var.type = "Hopf"
    free_var.ω = Num(0)
    free_var.symbol = new_symbol
    free_var.name = var_name(new_symbol)

    # define Problem as usual but with the Hopf Jacobian (always computed implicitly)
    p = Problem(eom; Jacobian=J)
    return p
end

