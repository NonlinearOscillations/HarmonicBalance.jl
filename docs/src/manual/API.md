# [API](@id doc-API)

**Table of contents**

[[toc]] <!-- the level setting is in ".vitepress/config.mts" -->

## System objects and types

```@docs
HarmonicBalance.d
DifferentialEquation
HarmonicVariable
HarmonicEquation
```

```@docs
rearrange_standard
rearrange_standard!
first_order_transform!
is_rearranged_standard
get_equations
```

```@docs
get_harmonic_equations
get_krylov_equations
add_harmonic!
```

```@docs
get_independent_variables
get_variables
```

## Solving and transforming solutions

```@docs
get_steady_states
```

### Methods

```@docs
WarmUp
TotalDegree
Polyhedral
```

### Access solutions

```@docs
get_single_solution
transform_solutions
```

### Classify

```@docs
classify_solutions!
get_class
```

### Plotting

```@docs
plot
plot!
plot_phase_diagram
plot_spaghetti
```

## Limit cycles

```@docs
get_limit_cycles
get_cycle_variables
add_pairs!
```

## Linear Response

```@autodocs
Modules = [HarmonicBalance.LinearResponse]
Private = false
Order = [:function]
```

```@docs
get_Jacobian
```

## Extensions

### OrdinaryDiffEq

```@docs
AdiabaticSweep
follow_branch
plot_1D_solutions_branch
```

```@autodocs; canonical=false
Modules = [Base.get_extension(HarmonicBalance, :TimeEvolution)]
Private = false
Order = [:function]
```

### SteadyStateSweep

```@docs
steady_state_sweep
```

### ModelingToolkit

```@docs
ODEProblem
ModelingToolkit.ODESystem
ModelingToolkit.SciMLBase.SteadyStateProblem
ModelingToolkit.SciMLBase.NonlinearProblem
```
