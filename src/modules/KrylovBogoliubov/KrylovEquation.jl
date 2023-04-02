using HarmonicBalance: flatten

#TODO Can KB have variables with zero harmonics?
# function van_der_Pol(diff_eom::DifferentialEquation, time::Num)
#     !is_harmonic(diff_eom, time) && error("The differential equation is not harmonic in ", time, " !")
#     eqs = collect(values(diff_eom.equations))
#     rules, vars = Dict(), []

#     # keep count to label new variables
#     uv_idx = 1;

#     ω = values(diff_eom.harmonics) |> unique |> flatten |> first

#     for nvar in get_variables(diff_eom) # sum over natural variables
#         to_substitute = Num(0) # combine all the subtitution rules for var

#         rule_u, hvar_u = _create_harmonic_variable(nvar, ω, time, "u", new_symbol="u"*string(uv_idx))
#         rule_v, hvar_v = _create_harmonic_variable(nvar, ω, time, "v", new_symbol="v"*string(uv_idx))
#         rule = rule_u - rule_v
#         push!(vars, hvar_u, hvar_v)

#         to_substitute += rule


#         rules[nvar] = to_substitute # total sub rule for nvar
#     end
#     eqs = substitute_all(eqs, rules)
#     HarmonicEquation(eqs, Vector{HarmonicVariable}(vars), diff_eom)
# end



function get_krylov_equations(diff_eom::DifferentialEquation; fast_time=nothing, slow_time=nothing)

    slow_time = isnothing(slow_time) ? (@variables T; T) : slow_time
    fast_time = isnothing(fast_time) ? get_independent_variables(diff_eom)[1] : fast_time

    harmonics = values(sys.harmonics)
    all(isempty.(harmonics)) && error("No harmonics specified!")
    any(isempty.(harmonics)) && error("Krylov-Bogoliubov method needs all vairables to have a single harmonic!")
    any(length.(harmonics) .> 1) && error("Krylov-Bogoliubov method only supports a single harmonic!")

    # eom = harmonic_ansatz(diff_eom, fast_time); # substitute trig functions into the differential equation
    # eom = slow_flow(eom, fast_time=fast_time, slow_time=slow_time; degree=2); # drop 2nd order time derivatives
    # fourier_transform!(eom, fast_time); # perform averaging over the frequencies originally specified in dEOM
    # ft_eom_simplified = drop_powers(eom, d(get_variables(eom), slow_time), 2); # drop higher powers of the first-order derivatives
    # return ft_eom_simplified
end
