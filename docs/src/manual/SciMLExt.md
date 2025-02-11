# Extention to the SciML ecosystem

The [SciML ecosystem](https://sciml.ai/) provides a rich set of tools to solve (non)-linear equations, differential equations and inverse problems. We provide an interface (in the form of [Package extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions))) to export the the derived harmonic equations computed with [harmonic balance method](@ref Harmonic_Balance) or [krylov-bogoliubov method](@ref Krylov-Bogoliubov) to the SciML ecosystem.

## ModeligToolkit.jl

The [`ModelingToolkit.jl`](https://github.com/SciML/ModelingToolkit.jl) (MTK) package provides a symbolic framework for defining and simplifying mathematical models. Through, MTK SciML provides a symbolic interface for their ecosystem

```@docs; canonical=false
ODEProblem(eom::Union{DifferentialEquation, HarmonicEquation},
    u0,
    tspan::Tuple,
    p::AbstractDict;
    kwargs...
)
ModelingToolkit.ODESystem
ModelingToolkit.SciMLBase.SteadyStateProblem
ModelingToolkit.SciMLBase.NonlinearProblem
```
