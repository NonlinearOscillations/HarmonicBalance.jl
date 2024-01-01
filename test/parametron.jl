using HarmonicBalance
using Symbolics
#using Test # do not use Test as this file is used for precompilation

@variables Ω γ λ F x θ η α ω0 ω t T ψ
@variables x(t)

natural_equation = d(d(x,t),t) + γ*d(x,t) + Ω^2*(1-λ*cos(2*ω*t+ψ))*x + α*x^3 + η*d(x,t)*x^2
forces = F*cos(ω*t+θ)
dEOM = DifferentialEquation(natural_equation + forces, x)
add_harmonic!(dEOM, x, ω)
harmonic_eq = get_harmonic_equations(dEOM, slow_time=T, fast_time=t);
p = HarmonicBalance.Problem(harmonic_eq);


fixed = (Ω => 1.0,γ => 1E-2, λ => 5E-2, F => 1E-3,  α => 1.,  η=>0.3, θ => 0, ψ => 0)
varied = ω => LinRange(0.9, 1.1, 100)
res = get_steady_states(p, varied, fixed, show_progress=false);

p = HarmonicBalance.Problem(harmonic_eq, Jacobian="implicit");
res = get_steady_states(p, varied, fixed, show_progress=false);

classify_solutions!(res, "sqrt(u1^2 + v1^2) > 1E-10", "nonzero")

# save the result, try and load in the next step
# current_path = @__DIR__
# HarmonicBalance.save(current_path * "/parametron_result.jld2", res)

# try to run a 2D calculation
fixed = (Ω => 1.0,γ => 1e-2, F => 1e-3, α => 1.0, η=>0.3, θ => 0, ψ => 0)
varied = (ω => LinRange(0.9, 1.1, 10), λ => LinRange(0.01, 0.05, 10))
res = get_steady_states(p, varied, fixed, show_progress=false);


###
# COMPARE TO KNOWN RESULTS
###
@variables u1, v1
ref1 = (Ω^2)*u1 + F*cos(θ) + γ*Differential(T)(u1) + (3//4)*α*(u1^3) + γ*ω*v1 + (2//1)*ω*Differential(T)(v1) + (1//4)*η*ω*(v1^3) + (3//4)*η*(u1^2)*Differential(T)(u1) + (1//4)*η*(v1^2)*Differential(T)(u1) + (3//4)*α*(v1^2)*u1 + (1//4)*η*ω*(u1^2)*v1 + (1//2)*η*u1*v1*Differential(T)(v1) + (1//2)*λ*(Ω^2)*v1*sin(ψ) - (ω^2)*u1 - (1//2)*λ*(Ω^2)*u1*cos(ψ)
ref2 = γ*Differential(T)(v1) + (Ω^2)*v1 + (3//4)*α*(v1^3) + (3//4)*α*(u1^2)*v1 + (1//4)*η*(u1^2)*Differential(T)(v1) + (3//4)*η*(v1^2)*Differential(T)(v1) + (1//2)*λ*(Ω^2)*v1*cos(ψ) + (1//2)*η*u1*v1*Differential(T)(u1) + (1//2)*λ*(Ω^2)*u1*sin(ψ) - F*sin(θ) - (ω^2)*v1 - (2//1)*ω*Differential(T)(u1) - γ*ω*u1 - (1//4)*η*ω*(u1^3) - (1//4)*η*ω*(v1^2)*u1

averaged = HarmonicBalance._remove_brackets(harmonic_eq)
@assert isequal(simplify(expand(averaged[1]-ref1)), 0)
@assert isequal(simplify(expand(averaged[2]-ref2)), 0)
