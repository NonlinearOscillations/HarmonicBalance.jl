![tests](https://github.com/NonlinearOscillations/HarmonicBalance.jl/workflows/Run%20tests/badge.svg?branch=master)

# HarmonicBalance.jl

**HarmonicBalance.jl** is a Julia package for solving nonlinear differential equations using the method of harmonic balance.

## Installation

To install HarmonicBalance.jl, you can use the github repo or the Julia package manager,
```julia
using Pkg
Pkg.add("HarmonicBalance")
```

## Documentation

For a detailed description of the package and examples, see the [stable documentation](https://arxiv.org/).

[This repo](https://github.com/NonlinearOscillations/HarmonicBalance-notebooks) contains a collection of example notebooks.

## Simple example
Let's find the steady states of a driven Duffing oscillator with nonlinear damping, its equation of motion is:

<img src="/docs/images/DuffingEq.png" width="450">

```julia
using HarmonicBalance
@variables α, ω, ω0, F, t, T, η, x(t) # declare constant variables and a function x(t)
diff_eq = DifferentialEquation(d(x,t,2) + ω0*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t), x)
add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq, slow_time=T, fast_time=t)

problem = Problem(harmonic_eq) # a steady-state problem

fixed = ParameterList(α => 1., ω0 => 1.0, F => 0.01, η=>0.1)   # fixed parameters
swept = ParameterRange(ω => LinRange(0.9, 1.2, 100))         # range of parameter values
solutions = get_steady_states(problem, swept, fixed)
```
```
A steady state result for 100 parameter points

Solution branches:   3
   of which real:    3
   of which stable:  2

Classes: stable, physical, Hopf, binary_labels
```

```julia
plot_1D_solutions(solutions, x="ω", y="sqrt(u1^2 + v1^2)"]);
```

<img src="/docs/images/DuffingPlot.png" width="500">

## Citation

If you use HarmonicBalance.jl in your project, we kindly ask you to cite [this paper](https://arxiv.org/):

```bib
@article{}
}
```

