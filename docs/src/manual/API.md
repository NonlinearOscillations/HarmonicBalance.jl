```@meta
CurrentModule = HarmonicBalance
```

# [API](@id doc-API)

**Table of contents**

[[toc]] <!-- the level setting is in ".vitepress/config.mts" -->

## System objects and types

```@docs
DifferentialEquation
HarmonicVariable
HarmonicEquation
```
```@docs; canonical=false
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

```@docs
get_Jacobian
```

```@docs
plot_eigenvalues
plot_linear_response
plot_rotframe_jacobian_response
```

## Extentsions

### OrdinaryDiffEq
```@docs
AdiabaticSweep
follow_branch
plot_1D_solutions_branch
```

### SteadyStateSweep
```@docs
steady_state_sweep
```

### ModelingToolkit
```@docs
```
