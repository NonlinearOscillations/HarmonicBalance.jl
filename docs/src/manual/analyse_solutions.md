# Extracting Solutions

After computing the steady-states of the harmonic equations, you'll want to extract the solutions from the [`HarmonicBalance.Result`](@ref) struct.

## Basic Solution Extraction

For plotting, you can extract the solutions using the `get_solutions` function, which parses a string into a symbolic expression, evaluates it for every steady state solution and filters the solutions by the requsted class.

```@docs; canonical=false
get_solutions
get_single_solution
```

## Attractors

```@docs; canonical=false
attractors
phase_diagram
```
