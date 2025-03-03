
# Methods

We offer several methods for solving the nonlinear algebraic equations that arise from the harmonic balance procedure. Each method has different tradeoffs between speed, robustness, and completeness.

## Total Degree Method {#Total-Degree-Method}
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.TotalDegree-manual-methods' href='#HarmonicBalance.TotalDegree-manual-methods'><span class="jlbinding">HarmonicBalance.TotalDegree</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TotalDegree
```


The Total Degree homotopy method performs a homotopy $H(x, t) = γ t G(x) + (1-t) F(x)$ from the trivial polynomial system $F(x) =xᵢ^{dᵢ} +aᵢ$ with the maximal degree $dᵢ$ determined by the [Bezout bound](https://en.wikipedia.org/wiki/B%C3%A9zout%27s_theorem). The method guarantees to find all solutions, however, it comes with a high computational cost. See [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/guides/totaldegree/) for more information.

**Fields**
- `gamma::Complex`: Complex multiplying factor of the start system G(x) for the homotopy
  
- `thread::Bool`: Boolean indicating if threading is enabled.
  
- `tracker_options::HomotopyContinuation.TrackerOptions`: Options for the tracker.
  
- `endgame_options::HomotopyContinuation.EndgameOptions`: Options for the endgame.
  
- `compile::Union{Bool, Symbol}`: Compilation options.
  
- `seed::UInt32`: Seed for random number generation.
  


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/methods.jl#L8-L20)

</details>


## Polyhedral Method {#Polyhedral-Method}
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.Polyhedral-manual-methods' href='#HarmonicBalance.Polyhedral-manual-methods'><span class="jlbinding">HarmonicBalance.Polyhedral</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
Polyhedral
```


The Polyhedral homotopy method constructs a homotopy based on the polyhedral structure of the polynomial system. It is more efficient than the Total Degree method for sparse systems, meaning most of the coefficients are zero. It can be especially useful if you don&#39;t need to find the zero solutions (`only_non_zero = true`), resulting in a speed up. See [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/guides/polyhedral/) for more information.

**Fields**
- `only_non_zero::Bool`: Boolean indicating if only non-zero solutions are considered.
  
- `thread::Bool`: Boolean indicating if threading is enabled.
  
- `tracker_options::HomotopyContinuation.TrackerOptions`: Options for the tracker.
  
- `endgame_options::HomotopyContinuation.EndgameOptions`: Options for the endgame.
  
- `compile::Union{Bool, Symbol}`: Compilation options.
  
- `seed::UInt32`: Seed for random number generation.
  


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/methods.jl#L57-L69)

</details>


## Warm Up Method {#Warm-Up-Method}
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.WarmUp-manual-methods' href='#HarmonicBalance.WarmUp-manual-methods'><span class="jlbinding">HarmonicBalance.WarmUp</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
WarmUp
```


The Warm Up method prepares a warmup system with the Total Degree method using the parameter at `index` perturbed by `perturbation_size`. The warmup system is used to perform a homotopy using all other systems in the parameter sweep. It is very efficient for systems with minimal bifurcation in the parameter sweep. The Warm Up method should in theory guarantee to find all solutions, however, if the `start_parameters` is not proper (to close to the real line) it could miss some solutions.

See[HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/guides/many-systems/) for more information.

**Fields**
- `warm_up_method::Union{Polyhedral{T}, TotalDegree{T}} where T`: Method used for the warmup system.
  
- `start_parameters::Vector`: Start parameters.
  
- `thread::Bool`: Boolean indicating if threading is enabled.
  
- `tracker_options::HomotopyContinuation.TrackerOptions`: Options for the tracker.
  
- `endgame_options::HomotopyContinuation.EndgameOptions`: Options for the endgame.
  
- `compile::Union{Bool, Symbol}`: Compilation options.
  
- `seed::UInt32`: Seed for random number generation.
  


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/methods.jl#L110-L125)

</details>

