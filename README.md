<img src="./logo.png" width="750" align="top">

![tests](https://github.com/NonlinearOscillations/HarmonicBalance.jl/workflows/Run%20tests/badge.svg?branch=master)
[![docs](https://img.shields.io/badge/docs-online-blue.svg)](https://nonlinearoscillations.github.io/HarmonicBalance.jl/stable/)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FHarmonicBalance&query=total_requests&label=Downloads)](https://juliapkgstats.com/pkg/HarmonicBalance)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FHarmonicBalance&query=total_requests&suffix=%2Fmonth&label=Downloads)](https://juliapkgstats.com/pkg/HarmonicBalance)

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%E2%9C%88%EF%B8%8F%20tested%20with%20-%20JET.jl%20-%20red)](https://github.com/aviatesk/JET.jl)


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
@variables α ω ω0 F η t x(t) # declare constant variables and a function x(t)
diff_eq = DifferentialEquation(d(x,t,2) + ω0^2*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t), x)
add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq)

# solve for a range of ω
result = get_steady_states(harmonic_eq, (ω => range(0.9, 1.2, 100),
                                          α => 1., ω0 => 1.0, F => 0.01, η => 0.1))
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
