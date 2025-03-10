<img src="./logo.png" width="750" align="top">

[![docs](https://img.shields.io/badge/docs-online-blue.svg)](https://quantumengineeredsystems.github.io/HarmonicBalance.jl/)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FHarmonicBalance&query=total_requests&label=Downloads)](https://juliapkgstats.com/pkg/HarmonicBalance)
[![Downloads](https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Fmonthly_downloads%2FHarmonicBalance&query=total_requests&suffix=%2Fmonth&label=Downloads)](https://juliapkgstats.com/pkg/HarmonicBalance)

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%E2%9C%88%EF%B8%8F%20tested%20with%20-%20JET.jl%20-%20red)](https://github.com/aviatesk/JET.jl)


**HarmonicBalance.jl** is a Julia package for solving periodic, nonlinear differential equations using the method of harmonic balance.

## Installation

To install HarmonicBalance.jl, you can use the github repo or the Julia package manager,
```julia
using Pkg
Pkg.add("HarmonicBalance")
```

## Documentation

For a detailed description of the package and examples, see the [documentation](https://quantumengineeredsystems.github.io/HarmonicBalance.jl).

[This repo](https://github.com/quantumengineeredsystems/HarmonicBalance-notebooks) contains a collection of example notebooks.

## [Example: steady states of a nonlinear resonator](https://quantumengineeredsystems.github.io/HarmonicBalance.jl/dev/tutorials/steady_states)
Let's find the steady states of a driven Duffing resonator with nonlinear damping, with equation of motion:

<img src="/docs/images/github_readme_eq.png" width="450">

```julia
using HarmonicBalance, Plots
@variables α ω ω0 F η t x(t) # declare constant variables and a function x(t)
diff_eq = DifferentialEquation(d(x,t,2) + ω0^2*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t), x)
add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq)

# solve for a range of ω
result = get_steady_states(harmonic_eq, ω => range(0.9, 1.2, 100), (α => 1., ω0 => 1.0, F => 0.01, η => 0.1))
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

### BibTeX entry

```bibtex
@Article{10.21468/SciPostPhysCodeb.6,
    title={{HarmonicBalance.jl: A Julia suite for nonlinear dynamics using harmonic balance}},
    author={Jan Košata and Javier del Pino and Toni L. Heugel and Oded Zilberberg},
    journal={SciPost Phys. Codebases},
    pages={6},
    year={2022},
    publisher={SciPost},
    doi={10.21468/SciPostPhysCodeb.6},
    url={https://scipost.org/10.21468/SciPostPhysCodeb.6}
}
```

## See also

- [JosephsonCircuits.jl](https://github.com/kpobrien/JosephsonCircuits.jl): Models superconducting circuits using modified nodal analysis and harmonic balance with an analytic Jacobian.
- [Manlab](https://manlab.lma.cnrs-mrs.fr/spip/): A similar package in Matlab also using continuation methods and using the Harmonic Balance method for periodic orbits.
