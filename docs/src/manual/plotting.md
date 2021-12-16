# Analysis and plotting

The key method for visualization is `transform_solutions`, which parses a string into a symbolic expression and evaluates it for every steady state solution. 

```@docs
HarmonicBalance.transform_solutions
```

## Plotting solutions

Any function of the steady state solutions may be plotted. 
In 1D, the solutions are colour-coded according to the branches obtained by `sort_solutions`. 

```@docs
HarmonicBalance.plot_1D_solutions
HarmonicBalance.plot_1D_jacobian_eigenvalues
HarmonicBalance.plot_2D_solutions
```

## Plotting phase diagrams

In many problems, rather than in any property of the solutions themselves, we are interested in the phase diagrams, encoding the number of (stable) solutions in different regions of the parameter space.

```@docs
HarmonicBalance.plot_2D_phase_diagram
HarmonicBalance.plot_2D_phase_diagram_interactive
```