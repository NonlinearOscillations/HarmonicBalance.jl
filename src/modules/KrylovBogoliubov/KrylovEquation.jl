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
    eom.equations = average(eom, t)
end
function average(eom::HarmonicEquation, t)
    eqs = similar(eom.equations)
    for (i,eq) in pairs(eom.equations)
        lhs = average(Num(eq.lhs),t); rhs = average(Num(eq.rhs),t);
        eqs[i] = lhs ~ rhs
    end
    return eqs
end
function average(x, t)
    term = trig_reduce(x)
    indep = get_independent(term, t)
    ft = Num(simplify_complex(Symbolics.expand(indep)))
    Symbolics.expand(ft)
end

function reduce(x::Num)
    x = simplify_fractions.(x)
    x = trig_reduce.(Num.(x))
    x = Num.(simplify_complex.(expand.(x)))
end

function take_trig_integral(x::BasicSymbolic, ω, t)
    if isdiv(x)
        arg_num = arguments(x.num)
        return simplify(expand(sum(take_trig_integral.(arg_num, ω, t))*ω)) / (x.den*ω)
    else
        all_terms = get_all_terms(Num(x))
        trigs = filter(z -> is_trig(z), all_terms)
        D = Differential(t)

        rules = []
        for trig in trigs
            arg = arguments(trig.val) |> first
            type = operation(trig.val)

            term = (type == cos ? sin(arg) : -cos(arg)) / expand_derivatives(D(arg)) |> Num
            append!(rules, [trig => term])
        end

        return Symbolics.substitute(x, Dict(rules))
    end
end
take_trig_integral(x::Num, ω, t) =  take_trig_integral(Symbolics.expand(unwrap(x)), ω, t)


function get_krylov_equations(diff_eom::DifferentialEquation; order, fast_time=nothing, slow_time=nothing)

    order < 1  && error("The order of the Krylov-Bogoliubov method must be at least 1!")
    order > 2  && error("Krylov-Bogoliubov implemetation only supports up to second order!")

    slow_time = isnothing(slow_time) ? (@variables T; T) : slow_time
    fast_time = isnothing(fast_time) ? get_independent_variables(diff_eom)[1] : fast_time

    harmonics = values(diff_eom.harmonics)
    all(isempty.(harmonics)) && error("No harmonics specified!")
    any(isempty.(harmonics)) && error("Krylov-Bogoliubov method needs all vairables to have a single harmonic!")
    any(length.(harmonics) .> 1) && error("Krylov-Bogoliubov method only supports a single harmonic!")
    harmonics = first.(harmonics)

    deom = deepcopy(diff_eom)
    !is_rearranged_standard(deom) ? rearrange_standard!(deom) : nothing
    first_order_transform!(deom, fast_time)
    eom = van_der_Pol(deom, fast_time)

    eom = slow_flow(eom, fast_time=fast_time, slow_time=slow_time; degree=2)

    rearrange!(eom, d(get_variables(eom),T))
    eom.equations = expand.(simplify.(eom.equations))
    eom.equations = expand.(simplify.(eom.equations)) # need it two times to get it completely simplified

    if order == 1
        average!(eom, fast_time)
        return eom
    elseif order == 2
        vars_symb = get_variables(eom)
        Fₜ = Num.(getfield.(eom.equations, :lhs))
        F₀ = Num.(getfield.(average(eom, fast_time), :lhs))
        Fₜ′ = substitute(get_Jacobian(eom), Dict(zip(_remove_brackets.(vars_symb), vars_symb)))

        Ḋ₁ = reduce.(Fₜ - F₀)
        D₁ = take_trig_integral.(Ḋ₁, harmonics, fast_time)
        D₁ =  D₁ - average.(D₁, fast_time)

        Gₜ = reduce.(Fₜ′*D₁)
        G₀ = average.(Gₜ, fast_time)
        eom.equations = F₀ + G₀ .~ getfield.(eom.equations, :rhs)
    end

    return eom
end
