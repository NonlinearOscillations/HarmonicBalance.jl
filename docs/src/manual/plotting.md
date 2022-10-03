# Analysis and plotting

The key method for visualization is `transform_solutions`, which parses a string into a symbolic expression and evaluates it for every steady state solution. 

```@docs
HarmonicBalance.transform_solutions
```

## Plotting solutions

The function `plot` is multiple-dispatched to plot 1D and 2D datasets. 
In 1D, the solutions are colour-coded according to the branches obtained by `sort_solutions`. 

```@docs
HarmonicBalance.plot(::Result, varags...)
```



## Plotting phase diagrams

In many problems, rather than in any property of the solutions themselves, we are interested in the phase diagrams, encoding the number of (stable) solutions in different regions of the parameter space. `plot_phase_diagram` handles this for 1D and 2D datasets.

```@docs
HarmonicBalance.plot_phase_diagram
```