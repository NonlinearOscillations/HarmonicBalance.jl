<img src="https://user-images.githubusercontent.com/58295890/146517756-dc4692f7-6b79-4309-a04a-94008bab3517.png" width="750" align="top">

![tests](https://github.com/NonlinearOscillations/HarmonicBalance.jl/workflows/Run%20tests/badge.svg?branch=master)
[![docs](https://img.shields.io/badge/docs-online-blue.svg)](https://nonlinearoscillations.github.io/HarmonicBalance.jl/stable/)

**HarmonicBalance.jl** is a Julia package for solving nonlinear differential equations using the method of harmonic balance.

## Installation

To install HarmonicBalance.jl, you can use the github repo or the Julia package manager,
```julia
using Pkg
Pkg.add("HarmonicBalance")
```

## Documentation

For a detailed description of the package and examples, see the [documentation](https://nonlinearoscillations.github.io/HarmonicBalance.jl/).

[This repo](https://github.com/NonlinearOscillations/HarmonicBalance-notebooks) contains a collection of example notebooks.

## [Simple example](https://nonlinearoscillations.github.io/HarmonicBalance.jl/stable/examples/simple_Duffing/)
Let's find the steady states of a driven Duffing oscillator with nonlinear damping, its equation of motion is:

<img src="/docs/images/github_readme_eq.png" width="450">

```julia
using HarmonicBalance
@variables α, ω, ω0, F, t, η, x(t) # declare constant variables and a function x(t)
diff_eq = DifferentialEquation(d(x,t,2) + ω0^2*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t), x)
add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = (α => 1., ω0 => 1.0, F => 0.01, η=>0.1)   # fixed parameters
varied = ω => LinRange(0.9, 1.2, 100)           # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed)
```
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

<img src="/docs/images/github_readme_plot.png">

## Citation

If you use HarmonicBalance.jl in your project, we kindly ask you to cite [this paper](https://scipost.org/SciPostPhysCodeb.6), namely:

**HarmonicBalance.jl: A Julia suite for nonlinear dynamics using harmonic balance**
Jan Košata, Javier del Pino, Toni L. Heugel, Oded Zilberberg
SciPost Phys. Codebases 6 (2022) 
