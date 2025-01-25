function add_pairs!(eom::DifferentialEquation, ω_lc::Num)
    for var in get_variables(eom), ω in eom.harmonics[var]
        add_harmonic!(eom, var, ω + ω_lc)
        add_harmonic!(eom, var, ω - ω_lc)
    end
    return nothing
end

"""
    add_pairs!(eom::DifferentialEquation; ω_lc::Num, n=1)

Add a limit cycle harmonic `ω_lc` to the system
Equivalent to adding `n` pairs of harmonics ω +- ω_lc for each existing ω.
"""
function add_pairs!(eom::DifferentialEquation; ω_lc::Num, n::Int)
    foreach(1:n) do k
        add_pairs!(eom, ω_lc)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return the harmonic variables which participate in the limit cycle labelled by `ω_lc`.
"""
function get_cycle_variables(eom::HarmonicEquation, ω_lc::Num)
    return vars = filter(x -> any(isequal.(ω_lc, get_all_terms(x.ω))), eom.variables)
end

"""
    Obtain the Jacobian of `eom` with a gauge-fixed variable `fixed_var`.
    `fixed_var` marks the variable which is fixed to zero due to U(1) symmetry.
    This leaves behind a finite degeneracy of solutions (see Chapter 6 of Jan's thesis).

    For limit cycles, we always use an 'implicit' Jacobian - a function which only returns the numerical Jacobian when a numerical solution
    is inserted. Finding the analytical Jacobian is usually extremely time-consuming.
"""
function _gaugefixed_Jacobian(
    eom::HarmonicEquation, fixed_var::HarmonicVariable; sym_order, rules
)
    rules = Dict(rules)
    setindex!(rules, 0, _remove_brackets(fixed_var))
    jac = get_implicit_Jacobian(eom; rules=rules, sym_order=sym_order)
    return jac
end

"""
$(TYPEDSIGNATURES)

Construct a `Problem` from `eom` in the case where U(1) symmetry is present
due to having added a limit cycle frequency `ω_lc`.
"""
function _cycle_Problem(eom::HarmonicEquation, swept, fixed, ω_lc::Num)
    eom = deepcopy(eom) # do not mutate eom
    isempty(get_cycle_variables(eom, ω_lc)) ? error("No Hopf variables found!") : nothing
    if !any(isequal.(eom.parameters, ω_lc))
        error(ω_lc, " is not a parameter of the harmonic equation!")
    else
        nothing
    end

    # eliminate one of the u,v variables (gauge fixing)
    fixed_var = _choose_fixed(eom, ω_lc) # remove the HV of the lowest subharmonic

    # get the Hopf Jacobian before altering anything - this is the usual Jacobian but the entries corresponding
    # to the fixed variable are removed
    _fix_gauge!(eom, ω_lc, fixed_var)

    # define Problem as usual but with the Hopf Jacobian (always computed implicitly)
    p = HomotopyContinuationProblem(eom, swept, fixed; compute_Jacobian=false)
    return p
end

function _choose_fixed(eom, ω_lc)
    vars = get_cycle_variables(eom, ω_lc)
    return first(vars) # This is arbitrary; better would be to substitute with values
end

function limit_cycle_problem(eom::HarmonicEquation, swept::OrderedDict, fixed::OrderedDict, ω_lc)
    swept, fixed = promote_types(swept, fixed)
    prob = _cycle_Problem(eom, swept, fixed, ω_lc) # the eom in this problem is gauge-fixed
    jacobian = _gaugefixed_Jacobian(
        eom, # this is not gauge_fixed
        _choose_fixed(eom, ω_lc);
        sym_order=_free_symbols(prob),
        rules=fixed,
    )
    return HomotopyContinuationProblem(
        prob.variables,
        prob.parameters,
        prob.swept_parameters,
        prob.fixed_parameters,
        prob.system,
        jacobian,
        prob.eom,
    )
end
function limit_cycle_problem(eom::HarmonicEquation, swept, fixed, ω_lc)
    return limit_cycle_problem(eom, OrderedDict(swept), OrderedDict(fixed), ω_lc)
end

"""
    get_limit_cycles(
        eom::HarmonicEquation, method::HarmonicBalanceMethod, swept, fixed, ω_lc; kwargs...)

Variant of `get_steady_states` for a limit cycle problem characterised by a Hopf frequency (usually called ω_lc)

Solutions with ω_lc = 0 are labelled unphysical since this contradicts the assumption of distinct harmonic variables corresponding to distinct harmonics.
"""
function get_limit_cycles(eom::HarmonicEquation, swept, fixed, ω_lc; kwargs...)
    prob = limit_cycle_problem(eom, swept, fixed, ω_lc)
    return get_limit_cycles(prob, WarmUp(), swept, fixed, ω_lc; kwargs...)
end
function get_limit_cycles(
    eom::HarmonicEquation, method::HarmonicBalanceMethod, swept, fixed, ω_lc; kwargs...
)
    prob = limit_cycle_problem(eom, swept, fixed, ω_lc)
    return get_limit_cycles(prob, method, swept, fixed, ω_lc; kwargs...)
end
function get_limit_cycles(eom::HarmonicEquation, pairs::Dict, ω_lc; kwargs...)
    swept = filter(x -> length(x[2]) > 1, pairs)
    fixed = filter(x -> length(x[2]) == 1, pairs)
    return get_limit_cycles(eom, swept, fixed, ω_lc; kwargs...)
end
function get_limit_cycles(
    eom::HarmonicEquation, method::HarmonicBalanceMethod, pairs::Dict, ω_lc; kwargs...
)
    swept = filter(x -> length(x[2]) > 1, pairs)
    fixed = filter(x -> length(x[2]) == 1, pairs)
    return get_limit_cycles(eom, method, swept, fixed, ω_lc; kwargs...)
end
function get_limit_cycles(
    prob::HomotopyContinuationProblem,
    method::HarmonicBalanceMethod,
    swept,
    fixed,
    ω_lc;
    classify_default=true,
    kwargs...,
)
    result = get_steady_states(prob, method; classify_default=classify_default, kwargs...)
    classify_default ? _classify_limit_cycles!(result, ω_lc) : nothing
    return result
end

# if abs(ω_lc) < tol, set all classifications to false
# TOLERANCE HARDCODED FOR NOW
function _classify_limit_cycles!(res::Result, ω_lc::Num)
    ω_lc_idx = findfirst(x -> isequal(x, ω_lc), res.problem.variables)
    for idx in CartesianIndices(res.solutions),
        c in filter(x -> x != "binary_labels", keys(res.classes))

        res.classes[c][idx] .*= abs.(getindex.(res.solutions[idx], ω_lc_idx)) .> 1e-10
    end

    classify_unique!(res, ω_lc)

    unique_stable = find_branch_order(
        map(.*, res.classes["stable"], res.classes["unique_cycle"])
    )

    # branches which are unique but never stable
    unique_unstable = setdiff(
        find_branch_order(map(.*, res.classes["unique_cycle"], res.classes["physical"])),
        unique_stable,
    )
    return order_branches!(res, vcat(unique_stable, unique_unstable)) # shuffle to have relevant branches first
end

"""
$(TYPEDSIGNATURES)

Fix the gauge in `eom` where `ω_lc` is the limit cycle frequency by constraining `fixed_var` to zero and promoting `ω_lc` to a variable.
"""
function _fix_gauge!(eom::HarmonicEquation, ω_lc::Num, fixed_var::HarmonicVariable)
    new_symbol = HarmonicBalance.declare_variable(
        string(ω_lc), first(get_independent_variables(eom))
    )
    rules = Dict(ω_lc => new_symbol, fixed_var.symbol => Num(0))
    eom.equations = expand_derivatives.(substitute_all(eom.equations, rules))
    eom.parameters = setdiff(eom.parameters, [ω_lc]) # ω_lc is now NOT a parameter anymore

    fixed_var.type = "Hopf"
    fixed_var.ω = Num(0) # not associated with any harmonic
    fixed_var.symbol = new_symbol
    return fixed_var.name = var_name(new_symbol)
end

_fix_gauge!(x::Num, fixed_var) = substitute_all(x, _remove_brackets(fixed_var) => 0)
