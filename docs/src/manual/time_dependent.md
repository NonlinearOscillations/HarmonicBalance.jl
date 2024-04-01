# Time evolution

Generally, solving the ODE of oscillatory systems in time requires numerically tracking the oscillations. This is a computationally expensive process; however, using the harmonic ansatz removes the oscillatory time-dependence. Simulating instead the harmonic variables of a `HarmonicEquation` is vastly more efficient - a steady state of the system appears as a fixed point in multidimensional space rather than an oscillatory function.

The module `TimeEvolution` is used to interface `HarmonicEquation` with the powerful solvers contained in `DifferentialEquations.jl`. Time-dependent parameter sweeps are defined using the object `ParameterSweep`.

```@docs
ODEProblem(::HarmonicEquation, ::Any; timespan::Tuple)
ParameterSweep(::Dict, ::Tuple)
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