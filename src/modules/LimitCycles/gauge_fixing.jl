export get_cycle_variables
export add_pairs!

using HarmonicBalance: is_rearranged, rearrange_standard, _remove_brackets
using HarmonicBalance.LinearResponse: get_implicit_Jacobian, get_Jacobian
import HarmonicBalance: is_stable, is_physical, is_Hopf_unstable, order_branches!, classify_binaries!, find_branch_order


function add_pairs!(eom::DifferentialEquation, Δω::Num)
    for var in get_variables(eom), ω in eom.harmonics[var]
        add_harmonic!(eom, var, ω + Δω)
        add_harmonic!(eom, var, ω - Δω)
    end
end


"""
    add_pairs!(eom::DifferentialEquation; Δω::Num, n=1)

Add a limit cycle harmonic `Δω` to the system
Equivalent to adding `n` pairs of harmonics ω +- Δω for each existing ω.
"""
add_pairs!(eom::DifferentialEquation; Δω::Num, n::Int) = [add_pairs!(eom, Δω) for k in 1:n]


"""
$(TYPEDSIGNATURES)

Return the harmonic variables which participate in the limit cycle labelled by `Δω`.
"""
get_cycle_variables(eom::HarmonicEquation, Δω::Num) = filter(x -> any(isequal.(Δω, get_all_terms(x.ω))), eom.variables)


"""
    Obtain the Jacobian of `eom` with a gauge-fixed variable `fixed_var`.
    `fixed_var` marks the variable which is fixed to zero due to U(1) symmetry.
    This leaves behind a finite degeneracy of solutions (see Chapter 6 of Jan's thesis).

    For limit cycles, we always use an 'implicit' Jacobian - a function which only returns the numerical Jacobian when a numerical solution
    is inserted. Finding the analytical Jacobian is usually extremely time-consuming.
"""
function _cycle_Jacobian(eom::HarmonicEquation, fixed_var::HarmonicVariable)
        J = get_implicit_Jacobian(eom, _remove_brackets(fixed_var)=>0)
end


"""
$(TYPEDSIGNATURES)

Construct a `Problem` from `eom` in the case where U(1) symmetry is present
due to having added a limit cycle frequency `Δω`.
`explicit_Jacobian=true` attempts to derive a symbolic Jacobian (usually not possible).
"""
function _cycle_Problem(eom::HarmonicEquation, Δω::Num; explicit_Jacobian=false)

    eom = deepcopy(eom) # do not mutate eom
    isempty(get_cycle_variables(eom, Δω)) ? error("No Hopf variables found!") : nothing
    !any(isequal.(eom.parameters, Δω)) ? error(Δω, " is not a parameter of the harmonic equation!") : nothing

    # eliminate one of the u,v variables (gauge fixing)
    fixed_var = get_cycle_variables(eom, Δω)[1]

    # get the Hopf Jacobian before altering anything - this is the usual Jacobian but the entries corresponding
    # to the fixed variable are removed
    J = explicit_Jacobian ? substitute_all(get_Jacobian(eom), _remove_brackets(fixed_var) => 0) : _cycle_Jacobian(eom, fixed_var)
    _fix_gauge!(eom, Δω, fixed_var)

    # define Problem as usual but with the Hopf Jacobian (always computed implicitly)
    p = Problem(eom; Jacobian=J)
    return p
end


"""
    get_limit_cycles(eom::HarmonicEquation, swept, fixed, Δω; kwargs...)

Variant of `get_steady_states` for a limit cycle problem characterised by a Hopf frequency (usually called Δω)

Solutions with Δω = 0 are labelled unphysical since this contradicts the assumption of distinct harmonic variables corresponding to distinct harmonics.
"""
function get_limit_cycles(eom::HarmonicEquation, swept, fixed, Δω; explicit_Jacobian=false, kwargs...)
    prob = _cycle_Problem(eom, Δω, explicit_Jacobian=explicit_Jacobian);
    result = get_steady_states(prob, swept, fixed; random_warmup=true, threading=true, classify_default=true, kwargs...)

    _classify_limit_cycles!(result, Δω)
    explicit_Jacobian || (result.jacobian = result.problem.jacobian)

    result
end

get_limit_cycles(eom::HarmonicEquation,swept,fixed; cycle_harmonic, kwargs...) = get_limit_cycles(eom, swept, fixed, cycle_harmonic; kwargs...)

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
    fixed_var.ω = Num(0) # not associated with any harmonic
    fixed_var.symbol = new_symbol
    fixed_var.name = var_name(new_symbol)
end
