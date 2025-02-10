# [Plotting solutions](@id plotting)

HarmonicBalance.jl comes with a plotting module `PlotsExt` that allows you to visualize the steady states in the [`HarmonicBalance.Result`](@ref). The module is conditionally loaded based on the `Plots.jl` package being loaded.

The function `plot` is multiple-dispatched to plot 1D and 2D datasets.
In 1D, the solutions are colour-coded according to the branches obtained by `sort_solutions`.

```@docs; canonical=false
plot(::HarmonicBalance.Result, varags...)
plot!
```

## Plotting phase diagrams

In many problems, rather than in any property of the solutions themselves, we are interested in the phase diagrams, encoding the number of (stable) solutions in different regions of the parameter space. `plot_phase_diagram` handles this for 1D and 2D datasets.

```@docs; canonical=false
HarmonicBalance.plot_phase_diagram
```

## Plot spaghetti plot

Sometimes, it is useful to plot the quadratures of the steady states (u, v) in function of a swept parameter. This is done with `plot_spaghetti`.

```@docs; canonical=false
HarmonicBalance.plot_spaghetti
```
