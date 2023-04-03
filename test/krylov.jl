using HarmonicBalance, Symbolics

@variables t T x(t) y(t) # symbolic variables
@variables ω ω0 γ F α λ ψ θ η

eq = [d(d(x,t),t) + γ*d(x,t) + ω0^2*(1-λ*cos(2*ω*t+ψ))*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t+θ)]

diff_eom = DifferentialEquation(eq, [x])

add_harmonic!(diff_eom, x, ω) # x will rotate at ω

harmonic_eq = get_krylov_equations(diff_eom)

fixed = (ω0 => 1.0, γ => 0.005, α => 1.0, η => 0, F=> 0.0, ψ => 0.0, θ => 0.0)
varied = (ω => range(0.99, 1.01, 100), λ => range(1e-6, 0.05, 100))
res = get_steady_states(harmonic_eq, varied, fixed, threading =true)

# plot 2D result
plot_phase_diagram(res)


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

get_Jacobian(harmonic_eq)

harmonic_eq.variables


Symbolics.jacobian(harmonic_eq.equations, get_variables(harmonic_eq))
