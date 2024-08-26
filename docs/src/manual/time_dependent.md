# Time evolution

Generally, solving the ODE of oscillatory systems in time requires numerically tracking the oscillations. This is a computationally expensive process; however, using the harmonic ansatz removes the oscillatory time-dependence. Simulating instead the harmonic variables of a `HarmonicEquation` is vastly more efficient - a steady state of the system appears as a fixed point in multidimensional space rather than an oscillatory function.

The Extention `TimeEvolution` is used to interface `HarmonicEquation` with the solvers contained in `OrdinaryDiffEq.jl`. Time-dependent parameter sweeps are defined using the object `ParameterSweep`. To use the `TimeEvolution` extension, one must first load the `OrdinaryDiffEq.jl` package.
```@docs
ODEProblem(::HarmonicEquation, ::Any; timespan::Tuple)
ParameterSweep
```

## Plotting

```@docs
HarmonicBalance.plot(::OrdinaryDiffEq.ODESolution, ::Any, ::HarmonicEquation)
```

## Miscellaneous
Using a time-dependent simulation can verify solution stability in cases where the Jacobian is too expensive to compute.

```@docs
HarmonicBalance.is_stable
```