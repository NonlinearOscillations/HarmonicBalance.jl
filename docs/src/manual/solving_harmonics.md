# Solving harmonic equations

Once an equation of motion has been defined and converted to a `HarmonicEquation`, we may use the homotopy continuation method (as implemented in HomotopyContinuation.jl) to find steady states. This means that,
having called `get_harmonic_equations`, we need to set all time-derivatives to zero and parse the resulting algebraic equations into a `Problem`.

`Problem` holds the steady-state equations, as well as the symbolic Jacobian which is needed for stability / linear response calculations. 

Once defined, a `Problem` can be solved for a set of input parameters using `get_steady_states` to obtain `Result`.

```@docs
Problem
get_steady_states
Result
```


## Classifying solutions
The solutions in `Result` are accompanied by similarly-sized boolean arrays stored in the dictionary `Result.classes`. The classes can be used by the plotting functions to show/hide/label certain solutions.

By default, classes "physical", "stable" and "binary\_labels" are created. User-defined classification is possible with `classify_solutions!`.

```@docs
HarmonicBalance.classify_solutions!
```

## Saving and Loading
Simulation data stored in `Result`, in addition to parameter values of the simulation, is stored for further use via `save_result`. This function is a simple wrapper for `DrWatson.jl`'s `safesave` function which prevents unwanted data overwriting.
```@docs
HarmonicBalance.save_result
HarmonicBalance.load_result
``` 