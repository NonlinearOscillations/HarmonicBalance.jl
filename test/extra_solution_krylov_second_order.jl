using DrWatson, Plots
using HarmonicBalance, Symbolics; HB = HarmonicBalance
using Symbolics.Rewriters
using Symbolics: BasicSymbolic, unwrap, arguments, isdiv, isadd, ismul

@variables t T x(t) y(t) # symbolic variables
@variables ω ω0 γ F α λ ψ θ η

eq = [d(d(x,t),t) + γ*d(x,t) + ω0^2*(1-λ*cos(2*ω*t))*x + α*x^3 #= + η*d(x,t)*x^2 ~ F*cos(ω*t+θ)=#]

diff_eom = DifferentialEquation(eq, [x])

add_harmonic!(diff_eom, x, ω) # x will rotate at ω

harmonic_eq1 = get_krylov_equations(diff_eom, order = 1)
harmonic_eq2 = get_krylov_equations(diff_eom, order = 2)

fixed = (ω0 => 1.0, γ => 0.005, α => 1.0, η => 0, F => 0.0, ψ => 0.0, θ => 0.0, λ => 0.03)
varied = (ω => range(0.80, 1.5, 100))
res1D1 = get_steady_states(harmonic_eq1, varied, fixed, threading =true, random_warmup=false)
res1D2 = get_steady_states(harmonic_eq2, varied, fixed, threading =true, random_warmup=false)

# plot 1D result
res1D1_plot = plot(res1D1, x = "ω", y = "√(u1^2+v1^2)", class="physical", not_class="stable")
res1D2_plot = plot(res1D2, x = "ω", y = "√(u1^2+v1^2)", class="stable")

fixed = (ω0 => 1.0, γ => 0.005, α => 1.0, η => 0, F => 0.0, ψ => 0.0, θ => 0.0)
varied = (ω => range(0.80, 1.05, 100), λ => range(1e-6, 0.2, 100))
res1 = get_steady_states(harmonic_eq1, varied, fixed, threading =true, random_warmup=false)
res2 = get_steady_states(harmonic_eq2, varied, fixed, threading =true, random_warmup=false)

# plot 2D result
p1 = plot_phase_diagram(res1, class="stable", title="First order KB")
p2 = plot_phase_diagram(res2, class="stable", title="Second order KB")

plot_title = "HarmonicBalance.jl: "*savename(Dict(ω0 => 1.0, γ => 0.005, α => 1.0))
p = plot(p1, p2, layout = (1,2), size = (1200,400), plot_title=plot_title, plot_titlevspan=0.1, plot_titlefontsize=15)
# savefig(p, "test/parametron_test.png")

res1D2.jacobian

norm.(getindex.(res1D2.solutions,8))
@show res1D2.classes["stable"][70][1:end];

eigvecs1 = []
eigvecs2 = []
for idx in 1:length(res1D2.solutions)
    λs1 = eigvals(real.(res1D2.jacobian(get_single_solution(res1D2, index=idx, branch=2))))
    λs2 = eigvals(real.(res1D2.jacobian(get_single_solution(res1D2, index=idx, branch=4))))
    push!(eigvecs1, λs1)
    push!(eigvecs2, λs2)
end

res1D2_plot = plot(res1D2, x = "ω", y = "√(u1^2+v1^2)", class="stable", legendposition=:best, title="Stable solutions: "*savename(Dict(ω0 => 1.0, γ => 0.005, α => 1.0, λ => 0.03)))
bif_plot_KB2 = plot_phase_diagram(res2, class="stable", title="Second order KB: "*savename(Dict(ω0 => 1.0, γ => 0.005, α => 1.0, λ => 0.03)))
hline!([0.03], color=:red, linestyle=:dash, label="λ = 0.03", legendposition=:best, linewidth=2)

freq = res1D2.swept_parameters[ω]
eig_plot1 = plot(freq, first.(imag.(eigvecs1)), color=1, xlabel="ω", label="Im(λ)", title="eigenvalues branch 3")
plot!(freq, last.(imag.(eigvecs1)), color=1, label=nothing)
plot!(freq, first.(real.(eigvecs1)), color=2, xlabel="ω",  label="Re(λ)")
plot!(freq, last.(real.(eigvecs1)), color=2, label=nothing, ylims=(-0.07,0.07))

eig_plot2 = plot(freq, first.(imag.(eigvecs2)), color=3, xlabel="ω", label="Im(λ)", title="eigenvalues branch 5")
plot!(freq, last.(imag.(eigvecs2)), color=3, label=nothing)
plot!(freq, first.(real.(eigvecs2)), color=4, xlabel="ω",  label="Re(λ)")
plot!(freq, last.(real.(eigvecs2)), color=4, label=nothing, ylims=(-0.07,0.07))


p = plot(res1D2_plot, bif_plot_KB2, eig_plot1, eig_plot2 , layout=(2,2), size=(1200,1000))
savefig(p, "test/extra_solution_second_order_KB.png")
