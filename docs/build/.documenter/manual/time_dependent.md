
# Time evolution {#Time-evolution}

Generally, solving the ODE of oscillatory systems in time requires numerically tracking the oscillations. This is a computationally expensive process; however, using the harmonic ansatz removes the oscillatory time-dependence. Simulating instead the harmonic variables of a `HarmonicEquation` is vastly more efficient - a steady state of the system appears as a fixed point in multidimensional space rather than an oscillatory function.

The extension `TimeEvolution` is used to interface `HarmonicEquation` with the solvers contained in `OrdinaryDiffEq.jl`. Time-dependent parameter sweeps are defined using the object `AdiabaticSweep`. To use the `TimeEvolution` extension, one must first load the `OrdinaryDiffEq.jl` package.
<details class='jldocstring custom-block' open>
<summary><a id='SciMLBase.ODEProblem-Tuple{HarmonicEquation, Any}-manual-time_dependent' href='#SciMLBase.ODEProblem-Tuple{HarmonicEquation, Any}-manual-time_dependent'><span class="jlbinding">SciMLBase.ODEProblem</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ODEProblem(
    eom::HarmonicEquation,
    fixed_parameters;
    sweep,
    u0,
    timespan,
    perturb_initial,
    kwargs...
)

```


Creates an ODEProblem object used by OrdinaryDiffEqTsit5.jl from the equations in `eom` to simulate time-evolution within `timespan`. `fixed_parameters` must be a dictionary mapping parameters+variables to numbers (possible to use a solution index, e.g. solutions[x][y] for branch y of solution x). If `u0` is specified, it is used as an initial condition; otherwise the values from `fixed_parameters` are used.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/TimeEvolution/ODEProblem.jl#L3-L9)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.AdiabaticSweep-manual-time_dependent' href='#HarmonicBalance.AdiabaticSweep-manual-time_dependent'><span class="jlbinding">HarmonicBalance.AdiabaticSweep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Represents a sweep of one or more parameters of a `HarmonicEquation`. During a sweep, the selected parameters vary linearly over some timespan and are constant elsewhere.

Sweeps of different variables can be combined using `+`.

**Fields**
- `functions::Dict{Num, Function}`: Maps each swept parameter to a function.
  

**Examples**

```julia
# create a sweep of parameter a from 0 to 1 over time 0 -> 100
julia> @variables a,b;
julia> sweep = AdiabaticSweep(a => [0., 1.], (0, 100));
julia> sweep[a](50)
0.5
julia> sweep[a](200)
1.0

# do the same, varying two parameters simultaneously
julia> sweep = AdiabaticSweep([a => [0.,1.], b => [0., 1.]], (0,100))
```


Successive sweeps can be combined,

```julia
sweep1 = AdiabaticSweep(ω => [0.95, 1.0], (0, 2e4))
sweep2 = AdiabaticSweep(λ => [0.05, 0.01], (2e4, 4e4))
sweep = sweep1 + sweep2
```


multiple parameters can be swept simultaneously,

```julia
sweep = AdiabaticSweep([ω => [0.95;1.0], λ => [5e-2;1e-2]], (0, 2e4))
```


and custom sweep functions may be used.

```julia
ωfunc(t) = cos(t)
sweep = AdiabaticSweep(ω => ωfunc)
```



[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/types.jl#L10-L49)

</details>


## Plotting
<details class='jldocstring custom-block' open>
<summary><a id='RecipesBase.plot-Tuple{ODESolution, Any, HarmonicEquation}-manual-time_dependent' href='#RecipesBase.plot-Tuple{ODESolution, Any, HarmonicEquation}-manual-time_dependent'><span class="jlbinding">RecipesBase.plot</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot(soln::ODESolution, f::String, harm_eq::HarmonicEquation; kwargs...)
```


Plot a function `f` of a time-dependent solution `soln` of `harm_eq`.

**As a function of time**

```
plot(soln::ODESolution, f::String, harm_eq::HarmonicEquation; kwargs...)
```


`f` is parsed by Symbolics.jl

**parametric plots**

```
plot(soln::ODESolution, f::Vector{String}, harm_eq::HarmonicEquation; kwargs...)
```


Parametric plot of f[1] against f[2]

Also callable as plot!


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/TimeEvolution/plotting.jl#L4-L22)

</details>


## Miscellaneous

Using a time-dependent simulation can verify solution stability in cases where the Jacobian is too expensive to compute.
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.is_stable-manual-time_dependent' href='#HarmonicBalance.is_stable-manual-time_dependent'><span class="jlbinding">HarmonicBalance.is_stable</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
is_stable(
    soln::OrderedCollections.OrderedDict,
    eom::HarmonicEquation;
    timespan,
    tol,
    perturb_initial
)

```


Numerically investigate the stability of a solution `soln` of `eom` within `timespan`. The initial condition is displaced by `perturb_initial`.

Return `true` the solution evolves within `tol` of the initial value (interpreted as stable).


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/TimeEvolution/ODEProblem.jl#L61)



```julia
is_stable(
    soln::OrderedCollections.OrderedDict,
    res::HarmonicBalance.Result;
    kwargs...
) -> Any

```


Returns true if the solution `soln` of the Result `res` is stable. Stable solutions are real and have all Jacobian eigenvalues Re(λ) &lt;= 0. `im_tol` : an absolute threshold to distinguish real/complex numbers. `rel_tol`: Re(λ) considered &lt;=0 if real.(λ) &lt; rel_tol*abs(λmax)


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/classification.jl#L77)

</details>

