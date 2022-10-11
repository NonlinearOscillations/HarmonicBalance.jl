export add_Hopf!
export get_Hopf_variables
export replace_Hopf_variable!

using HarmonicBalance: is_rearranged, rearrange_standard, _remove_brackets
using HarmonicBalance.LinearResponse: get_implicit_Jacobian, get_Jacobian
import HarmonicBalance: is_stable, is_physical, is_Hopf_unstable, order_branches!, classify_binaries!, find_branch_order
#import HarmonicBalance.get_steady_states; export get_steady_states


function add_Hopf!(eom::DifferentialEquation, Δω::Num)
    for var in get_variables(eom), ω in eom.harmonics[var]
        add_harmonic!(eom, var, ω + Δω)
        add_harmonic!(eom, var, ω - Δω)
    end
end


"""
    add_Hopf!(eom::DifferentialEquation; Δω::Num, n=1)

Add a Hopf degree of freedom to the system by adding `n` pairs of harmonics ω +- Δω for each existing ω
"""
add_Hopf!(eom::DifferentialEquation; Δω::Num, n::Int) = [add_Hopf!(eom, Δω) for k in 1:n]


"""
$(TYPEDSIGNATURES)

Return the harmonic variables which participate in the limit cycle labelled by `Δω`.
"""
get_Hopf_variables(eom::HarmonicEquation, Δω::Num) = filter(x -> any(isequal.(Δω, get_all_terms(x.ω))), eom.variables)


"""
    Obtain the Jacobian of `eom` with a gauge-fixed (Hopf) variable `fixed_var`.
    `fixed_var` marks the variable which is fixed to zero due to U(1) symmetry. 
    This leaves behind a finite degeneracy of solutions (see Chapter 6 of thesis).

    The Jacobian is always 'implicit' - a function which only returns the numerical Jacobian when a numerical solution
    is inserted. Finding the analytical Jacobian is usually extremely time-consuming.
"""
function _Hopf_implicit_Jacobian(eom::HarmonicEquation, fixed_var::HarmonicVariable)
        J = get_implicit_Jacobian(eom, _remove_brackets(fixed_var)=>0)
end


"""
$(TYPEDSIGNATURES)

Construct a `Problem` from `eom` in the case where U(1) symmetry is present 
due to having added a Hopf variable with frequency `Δω`.
"""
function _Hopf_Problem(eom::HarmonicEquation, Δω::Num; explicit_Jacobian=false)

    eom = deepcopy(eom) # do not mutate eom
    isempty(get_Hopf_variables(eom, Δω)) ? error("No Hopf variables found!") : nothing
    !any(isequal.(eom.parameters, Δω)) ? error(Δω, " is not a parameter of the harmonic equation!") : nothing 

    # eliminate one of the Cartesian variables
    fixed_var = get_Hopf_variables(eom, Δω)[1] 
    
    # get the Hopf Jacobian before altering anything - this is the usual Jacobian but the entries corresponding
    # to the free variable are removed
    J = explicit_Jacobian ? substitute_all(get_Jacobian(eom), _remove_brackets(fixed_var) => 0) : _Hopf_implicit_Jacobian(eom, fixed_var)
    _fix_gauge!(eom, Δω, fixed_var)
    
    # define Problem as usual but with the Hopf Jacobian (always computed implicitly)
    p = Problem(eom; Jacobian=J)
    return p
end


"""
    get_steady_states(eom::HarmonicEquation, swept, fixed, Δω; kwargs...)   

Variant of `get_steady_states` for a limit cycle problem characterised by a Hopf frequency (usually called Δω)

Solutions with Δω = 0 are labelled unphysical since this contradicts the assumption of distinct harmonic variables corresponding to distinct harmonics.
"""
function get_steady_states(eom::HarmonicEquation, swept, fixed, Δω; explicit_Jacobian=false, kwargs...)   
    prob = _Hopf_Problem(eom, Δω, explicit_Jacobian=explicit_Jacobian);
    result = HarmonicBalance.get_steady_states(prob, swept, fixed; random_warmup=true, threading=true, classify_default=true, kwargs...)

    _classify_limit_cycles!(result, Δω)

    result
end

get_steady_states(eom::HarmonicEquation,swept,fixed; Δω, kwargs...) = get_steady_states(eom, swept, fixed, Δω; kwargs...)

# if abs(Δω) < tol, set all classifications to false
# TOLERANCE HARDCODED FOR NOW
function _classify_limit_cycles!(res::Result, Δω::Num)
    Δω_idx = findfirst(x -> isequal(x, Δω), res.problem.variables)
    for idx in CartesianIndices(res.solutions), c in filter(x -> x != "binary_labels", keys(res.classes))
        res.classes[c][idx] .*= abs.(getindex.(res.solutions[idx], Δω_idx)) .> 1E-10
    end

    classify_unique!(res, Δω)

    unique_stable = find_branch_order(map(.*, res.classes["stable"], res.classes["unique"]))

    # branches which are unique but never stable
    unique_unstable = setdiff(find_branch_order(map(.*, res.classes["unique"], res.classes["physical"])), unique_stable)
    order_branches!(res, vcat(unique_stable, unique_unstable)) # shuffle to have relevant branches first
end


"""
$(TYPEDSIGNATURES)

Fix the gauge in `eom` where `Δω` is the limit cycle frequency by constraining `fixed_var` to zero and promoting `Δω` to a variable. 
"""
function _fix_gauge!(eom::HarmonicEquation, Δω::Num, fixed_var::HarmonicVariable)

    new_symbol = HarmonicBalance.declare_variable(string(Δω), first(get_independent_variables(eom)))
    rules = Dict(Δω => new_symbol, fixed_var.symbol => Num(0))
    eom.equations = expand_derivatives.(substitute_all(eom.equations, rules))
    eom.parameters = setdiff(eom.parameters, [Δω]) # Δω is now NOT a parameter anymore
    
    fixed_var.type = "Hopf"
    fixed_var.ω = Num(0)
    fixed_var.symbol = new_symbol
    fixed_var.name = var_name(new_symbol)
end


