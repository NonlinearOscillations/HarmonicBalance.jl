using HarmonicBalance, Symbolics; HB = HarmonicBalance
using Symbolics.Rewriters
using Symbolics: BasicSymbolic, unwrap, arguments, isdiv, isadd

@variables t T x(t) y(t) # symbolic variables
@variables ω ω0 γ F α λ ψ θ η

eq = [d(d(x,t),t) + γ*d(x,t) + ω0^2*x]

diff_eom = DifferentialEquation(eq, [x])

add_harmonic!(diff_eom, x, ω) # x will rotate at ω
rearrange_standard!(diff_eom)

first_order_transform!(diff_eom, t)
eom = van_der_Pol(diff_eom, t)

eom = slow_flow(eom, fast_time=t, slow_time=T; degree=2)

HB.rearrange!(eom, d(get_variables(eom),T))
eom.equations = expand.(simplify.(eom.equations))
eom.equations = expand.(simplify.(eom.equations))


function get_Jacobian(eom::HarmonicEquation)
    rearr = !HarmonicBalance.is_rearranged(eom) ? HarmonicBalance.rearrange_standard(eom) : eom
    lhs = _remove_brackets(rearr)
    vars = _remove_brackets.(eom.variables)

    get_Jacobian(lhs, vars)
end
function get_Jacobian(eqs::Vector{Num}, vars::Vector{Num})
    length(eqs) == length(vars) || error("Jacobians are only defined for square systems!")
    M = Matrix{Num}(undef, length(vars), length(vars))

    for idx in CartesianIndices(M)
        M[idx] = expand_derivatives(d(eqs[idx[1]], vars[idx[2]]))
    end
    M
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
    term = HB.trig_reduce(x)
    indep = HB.get_independent(term, t)
    ft = Num(HB.simplify_complex(Symbolics.expand(indep)))
    Symbolics.expand(ft)
end

expand_div(x::BasicSymbolic) = isdiv(x) ? sum(map(y-> y/x.den, arguments(x.num))) : x
expand_div(x::Num) = x |> unwrap |> expand_div |> Num
simplify_terms(x::BasicSymbolic) = isadd(x) ? sum(map(y-> simplify(y), arguments(x))) : x
simplify_terms(x::Num) = x |> unwrap |> simplify_terms |> Num

function take_trig_int(x::Num, t)
    x = expand_div(x)
    all_terms = get_all_terms(x)
    trigs = filter(z -> HB.is_trig(z), all_terms)
    D = Differential(t)

    rules = []
    for trig in trigs
        arg = arguments(trig.val) |> first
        type = operation(trig.val)

        if type == cos
            term = Num(sin(arg) / expand_derivatives(D(arg)))
        elseif type == sin
            term = Num(-cos(arg) / expand_derivatives(D(arg)))
        end
        append!(rules, [trig => term])
    end
    result = Symbolics.substitute(x, Dict(rules)) |> Num
    return result
end

Fₜ = getfield.(eom.equations, :lhs)
F₀ = getfield.(average(eom, t), :lhs)
Fₜ′ = get_Jacobian(eom)

Ḋ₁ = Fₜ - F₀
Ḋ₁ = simplify_fractions.(Ḋ₁)
Ḋ₁ = HB.trig_reduce.(Num.(Ḋ₁))
Ḋ₁ = Num.(HB.simplify_complex.(Symbolics.expand.(Ḋ₁)))
D₁ =  take_trig_int.(ft,t)
simplify_fractions(sum(arguments(D₁[1].val)[2:4]))
Gₜ = Fₜ′*D₁
G₀ = average.(Gₜ,t)
