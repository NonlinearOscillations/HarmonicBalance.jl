using HarmonicBalance: flatten, is_harmonic, _create_harmonic_variable, HarmonicEquation
using HarmonicBalance: trig_reduce, get_independent, simplify_complex
export van_der_Pol, get_krylov_equations, average!

#TODO Can KB have variables with zero harmonics?
function van_der_Pol(eom::DifferentialEquation, t::Num)
    !is_harmonic(eom, t) && error("The differential equation is not harmonic in ", t, " !")
    eqs = equations(eom)
    rules, vars = Dict(), []

    # keep count to label new variables
    uv_idx = 1;
    ω = values(eom.harmonics) |> unique |> flatten |> first
    nvars = get_variables(eom)
    nvars = nvars[length(nvars)÷2+1:end]

    for nvar in nvars # sum over natural variables
        rule_u, hvar_u = _create_harmonic_variable(nvar, ω, t, "u", new_symbol="u"*string(uv_idx))
        rule_v, hvar_v = _create_harmonic_variable(nvar, ω, t, "v", new_symbol="v"*string(uv_idx))
        rule = rule_u - rule_v; rules[nvar] = rule;

        D = Differential(t); nvar_t = diff2term(D(unwrap(nvar)));
        vdP_rules = Dict(D(hvar_u.symbol) => 0, D(hvar_v.symbol) => 0)
        rules[nvar_t] = substitute(expand_derivatives(D(rule)), vdP_rules)

        uv_idx += 1
        push!(vars, hvar_u, hvar_v)
    end
    eqs = expand_derivatives.(substitute_all(eqs, rules))
    HarmonicEquation(eqs, Vector{HarmonicVariable}(vars), eom)
end

function average!(eom::HarmonicEquation, t)
    eqs = similar(eom.equations)
    for (i,eq) in pairs(eom.equations)
        lhs = average(Num(eq.lhs),t); rhs = average(Num(eq.rhs),t);
        eqs[i] = lhs ~ rhs
    end
    eom.equations = eqs
end

function average(x, t)
    term = trig_reduce(x)
    indep = get_independent(term, t)
    ft = Num(simplify_complex(Symbolics.expand(indep)))
    Symbolics.expand(ft)
end

function get_krylov_equations(diff_eom::DifferentialEquation; fast_time=nothing, slow_time=nothing)

    slow_time = isnothing(slow_time) ? (@variables T; T) : slow_time
    fast_time = isnothing(fast_time) ? get_independent_variables(diff_eom)[1] : fast_time

    harmonics = values(diff_eom.harmonics)
    all(isempty.(harmonics)) && error("No harmonics specified!")
    any(isempty.(harmonics)) && error("Krylov-Bogoliubov method needs all vairables to have a single harmonic!")
    any(length.(harmonics) .> 1) && error("Krylov-Bogoliubov method only supports a single harmonic!")

    !is_rearranged_standard(diff_eom) ? rearrange_standard!(diff_eom) : nothing
    first_order_transform!(diff_eom, fast_time)
    eom = van_der_Pol(diff_eom, fast_time)

    eom = slow_flow(eom, fast_time=fast_time, slow_time=slow_time; degree=2)

    rearrange!(eom, d(get_variables(eom),T))
    eom.equations = expand.(simplify.(eom.equations))
    eom.equations = expand.(simplify.(eom.equations))

    average!(eom, fast_time)

    return eom
end
