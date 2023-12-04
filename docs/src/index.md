# HarmonicBalance.jl

HarmonicBalance.jl is a Julia package for solving nonlinear differential equations using the method of harmonic balance.

[This repo](https://github.com/NonlinearOscillations/HarmonicBalance-notebooks) contains a collection of example notebooks.

## Installation

To install HarmonicBalance.jl, you can use the github repo https://github.com/NonlinearOscillations/HarmonicBalance.jl or the Julia package manager,
```julia
using Pkg
Pkg.add("HarmonicBalance")
```

## Citation

If you use HarmonicBalance.jl in your project, we kindly ask you to cite [this paper](https://scipost.org/SciPostPhysCodeb.6):

```bib
@article{10.21468/SciPostPhysCodeb.6,
	title={{HarmonicBalance.jl: A Julia suite for nonlinear dynamics using harmonic  balance}},
	author={Jan Košata and Javier del Pino and Toni L. Heugel and Oded Zilberberg},
	journal={SciPost Phys. Codebases},
	pages={6},
	year={2022},
	doi={10.21468/SciPostPhysCodeb.6},
	url={https://scipost.org/10.21468/SciPostPhysCodeb.6},
}
```

## [Simple example](https://nonlinearoscillations.github.io/HarmonicBalance.jl/stable/examples/simple_Duffing/)
Let's find the steady states of a driven Duffing oscillator with nonlinear damping, its equation of motion is:
```math
\begin{equation} \label{eq:duffing}
\underbrace{\ddot{x}(t) + \gamma \dot{x}(t) + \omega_0^2 x(t)}_{\text{damped harmonic oscillator}} + \underbrace{\alpha x(t)^3}_{\text{Duffing coefficient}} = \underbrace{F \cos(\omega t)}_{\text{periodic drive}}
\end{equation}
```

```julia
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
The results are `show`n:
```
A steady state result for 100 parameter points

Solution branches:   3
   of which real:    3
   of which stable:  2

Classes: stable, physical, Hopf, binary_labels
```
```julia
plot(result, "sqrt(u1^2 + v1^2)")
```
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/simple_Duffing/response_single.png" alignment="center" \>
``` ⠀

## Documentation

This documentation was built using [Documenter.jl](https://github.com/JuliaDocs).

```@example
using Dates # hide
println("Documentation built $(Dates.now()) with Julia $(VERSION)") # hide
```