# Installation

It is easy to install HarmonicBalance.jl as we are registered in the Julia General registry.
You can simply run the following command in the Julia REPL:
```julia
julia> using Pkg
julia> Pkg.add("HarmonicBalance")
```
or 
```julia
julia> ] # `]` should be pressed
julia> Pkg.add("HarmonicBalance")
```
You can check which version you have installled with the command 
```julia
julia> ]
julia> status HarmonicBalance
```

# Getting Started

Let us find the steady states of an external driven Duffing oscillator with nonlinear damping. Its equation of motion is:
```math
\begin{equation} \label{eq:duffing}
\underbrace{\ddot{x}(t) + \gamma \dot{x}(t) + \omega_0^2 x(t)}_{\text{damped harmonic oscillator}} + \underbrace{\alpha x(t)^3}_{\text{Duffing coefficient}} = \underbrace{F \cos(\omega t)}_{\text{periodic drive}}
\end{equation}
```

```@example getting_started
using HarmonicBalance
@variables α ω ω0 F t η x(t) # declare constant variables and a function x(t)
diff_eq = DifferentialEquation(d(x,t,2) + ω0^2*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t), x)
add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = (α => 1.0, ω0 => 1.0, F => 0.01, η => 0.1)   # fixed parameters
varied = ω => range(0.9, 1.2, 100)           # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed)
```
The found steady states can be plotted as a function of the driving frequency:
```@example getting_started
plot(result, "sqrt(u1^2 + v1^2)")
```
