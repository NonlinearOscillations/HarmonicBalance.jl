# Solving harmonic equations

Once a differential equation of motion has been defined in `DifferentialEquation` and converted to a `HarmonicEquation`, we may use the homotopy continuation method (as implemented in HomotopyContinuation.jl) to find steady states. This means that,
having called `get_harmonic_equations`, we need to set all time-derivatives to zero and parse the resulting algebraic equations into a `Problem`.

`Problem` holds the steady-state equations, and (optionally) the symbolic Jacobian which is needed for stability / linear response calculations. 

Once defined, a `Problem` can be solved for a set of input parameters using `get_steady_states` to obtain `Result`.

```@docs; canonical=false
HarmonicBalance.Problem
get_steady_states
HarmonicBalance.Result
```


## Classifying solutions
The solutions in `Result` are accompanied by similarly-sized boolean arrays stored in the dictionary `Result.classes`. The classes can be used by the plotting functions to show/hide/label certain solutions.

By default, classes "physical", "stable" and "binary\_labels" are created. User-defined classification is possible with `classify_solutions!`.

```@docs
HarmonicBalance.classify_solutions!
```

## Sorting solutions
Solving a steady-state problem over a range of parameters returns a solution set for each parameter. For a continuous change of parameters, each solution in a set usually also changes continuously; it is said
to form a ''solution branch''. For an example, see the three colour-coded branches for the Duffing oscillator in Example 1.

For stable states, the branches describe a system's behaviour under adiabatic parameter changes. 

Therefore, after solving for a parameter range, we want to order each solution set such that the solutions' order reflects the branches.

The function `sort_solutions` goes over the the raw output of `get_steady_states` and sorts each entry such that neighboring solution sets minimize Euclidean distance.

Currently, `sort_solutions` is compatible with 1D and 2D arrays of solution sets.

```@docs
HarmonicBalance.sort_solutions
```