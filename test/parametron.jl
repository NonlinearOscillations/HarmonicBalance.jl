using HarmonicBalance
using Symbolics
using Test

@variables Ω,γ,λ,F, x,θ,η,α, ω0, ω,t,T, ψ
@variables x(t)

natural_equation =  d(d(x,t),t) + γ*d(x,t) + Ω^2*(1-λ*cos(2*ω*t+ψ))*x + α * x^3 +η *d(x,t) * x^2
forces =  F*cos(ω*t+θ)
dEOM = HarmonicBalance.DifferentialEquation(natural_equation + forces, x)
HarmonicBalance.add_harmonic!(dEOM, x, ω)
averagedEOM = HarmonicBalance.get_harmonic_equations(dEOM, slow_time=T, fast_time=t);
HarmonicBalance.Problem(averagedEOM);



###
# COMPARE TO KNOWN RESULTS
###
@variables u1, v1
ref1 = (Ω^2)*u1 + F*cos(θ) + γ*Differential(T)(u1) + (3//4)*α*(u1^3) + γ*ω*v1 + (2//1)*ω*Differential(T)(v1) + (1//4)*η*ω*(v1^3) + (3//4)*η*(u1^2)*Differential(T)(u1) + (1//4)*η*(v1^2)*Differential(T)(u1) + (3//4)*α*(v1^2)*u1 + (1//4)*η*ω*(u1^2)*v1 + (1//2)*η*u1*v1*Differential(T)(v1) + (1//2)*λ*(Ω^2)*v1*sin(ψ) - (ω^2)*u1 - (1//2)*λ*(Ω^2)*u1*cos(ψ)
ref2 = γ*Differential(T)(v1) + (Ω^2)*v1 + (3//4)*α*(v1^3) + (3//4)*α*(u1^2)*v1 + (1//4)*η*(u1^2)*Differential(T)(v1) + (3//4)*η*(v1^2)*Differential(T)(v1) + (1//2)*λ*(Ω^2)*v1*cos(ψ) + (1//2)*η*u1*v1*Differential(T)(u1) + (1//2)*λ*(Ω^2)*u1*sin(ψ) - F*sin(θ) - (ω^2)*v1 - (2//1)*ω*Differential(T)(u1) - γ*ω*u1 - (1//4)*η*ω*(u1^3) - (1//4)*η*ω*(v1^2)*u1

averaged = HarmonicBalance._remove_brackets(averagedEOM)
@test isequal(simplify(expand(averaged[1]-ref1)), 0)
@test isequal(simplify(expand(averaged[2]-ref2)), 0)
