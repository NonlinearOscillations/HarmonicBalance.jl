
# API {#doc-API}

**Table of contents**

[[toc]] &lt;!– the level setting is in &quot;.vitepress/config.mts&quot; –&gt;

## System objects and types {#System-objects-and-types}
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.d' href='#HarmonicBalance.d'><span class="jlbinding">HarmonicBalance.d</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



The derivative of f w.r.t. x of degree deg


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicVariable.jl#L108)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.DifferentialEquation' href='#HarmonicBalance.DifferentialEquation'><span class="jlbinding">HarmonicBalance.DifferentialEquation</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct DifferentialEquation
```


Holds differential equation(s) of motion and a set of harmonics to expand each variable. This is the primary input for `HarmonicBalance.jl`. After inputting the equations, the harmonics ansatz needs to be specified using `add_harmonic!`.

**Fields**
- `equations::OrderedCollections.OrderedDict{Num, Equation}`: Assigns to each variable an equation of motion.
  
- `harmonics::OrderedCollections.OrderedDict{Num, OrderedCollections.OrderedSet{Num}}`: Assigns to each variable a set of harmonics.
  

**Example**

```julia
julia> @variables t, x(t), y(t), ω0, ω, F, k;

# equivalent ways to enter the simple harmonic oscillator
julia> DifferentialEquation(d(x,t,2) + ω0^2 * x - F * cos(ω*t), x);
julia> DifferentialEquation(d(x,t,2) + ω0^2 * x ~ F * cos(ω*t), x);

# two coupled oscillators, one of them driven
julia> DifferentialEquation(
    [d(x,t,2) + ω0^2 * x - k*y, d(y,t,2) + ω0^2 * y - k*x] .~ [F * cos(ω*t), 0], [x,y]
);
```



[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/DifferentialEquation.jl#L1)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.HarmonicVariable' href='#HarmonicBalance.HarmonicVariable'><span class="jlbinding">HarmonicBalance.HarmonicVariable</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct HarmonicVariable
```


Holds a variable stored under `symbol` describing the harmonic `ω` of `natural_variable`.

**Fields**
- `symbol::Num`: Symbol of the variable in the HarmonicBalance namespace.
  
- `name::String`: Human-readable labels of the variable, used for plotting.
  
- `type::String`: Type of the variable (u or v for quadratures, a for a constant, Hopf for Hopf etc.)
  
- `ω::Num`: The harmonic being described.
  
- `natural_variable::Num`: The natural variable whose harmonic is being described.
  


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicVariable.jl#L1)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.HarmonicEquation' href='#HarmonicBalance.HarmonicEquation'><span class="jlbinding">HarmonicBalance.HarmonicEquation</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct HarmonicEquation
```


Holds a set of algebraic equations governing the harmonics of a `DifferentialEquation`.

**Fields**
- `equations::Vector{Equation}`: A set of equations governing the harmonics.
  
- `variables::Vector{HarmonicVariable}`: A set of variables describing the harmonics.
  
- `parameters::Vector{Num}`: The parameters of the equation set.
  
- `natural_equation::DifferentialEquation`: The natural equation (before the harmonic ansatz was used).
  


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicEquation.jl#L1)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.rearrange_standard' href='#HarmonicBalance.rearrange_standard'><span class="jlbinding">HarmonicBalance.rearrange_standard</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
rearrange_standard(
    eom::HarmonicEquation
) -> HarmonicEquation

```


Rearrange `eom` to the standard form, such that the derivatives of the variables are on one side.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicEquation.jl#L178)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.rearrange_standard!' href='#HarmonicBalance.rearrange_standard!'><span class="jlbinding">HarmonicBalance.rearrange_standard!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
rearrange_standard!(eom::DifferentialEquation)
rearrange_standard!(eom::DifferentialEquation, degree)

```


Rearranges the differential equations in `eom` to standard form, where the highest derivative of each variable (specified by `degree`, default 2) appears isolated on the left-hand side. Modifies the equations in place.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/DifferentialEquation.jl#L160)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.first_order_transform!' href='#HarmonicBalance.first_order_transform!'><span class="jlbinding">HarmonicBalance.first_order_transform!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
first_order_transform!(diff_eom::DifferentialEquation, time)

```


Transforms a higher-order differential equation system into an equivalent first-order system by introducing additional variables. Modifies the system in place. The `time` parameter specifies the independent variable used for differentiation.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/DifferentialEquation.jl#L202)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.is_rearranged_standard' href='#HarmonicBalance.is_rearranged_standard'><span class="jlbinding">HarmonicBalance.is_rearranged_standard</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
is_rearranged_standard(eom::DifferentialEquation) -> Any
is_rearranged_standard(
    eom::DifferentialEquation,
    degree
) -> Any

```


Checks if the differential equations in `eom` are arranged in standard form, where the highest derivative of each variable appears isolated on the left-hand side. The default degree is 2, corresponding to second-order differential equations.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/DifferentialEquation.jl#L147)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.get_equations' href='#HarmonicBalance.get_equations'><span class="jlbinding">HarmonicBalance.get_equations</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_equations(eom::DifferentialEquation) -> Vector{Equation}

```


Return the equations of `eom`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/DifferentialEquation.jl#L140)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.get_harmonic_equations' href='#HarmonicBalance.get_harmonic_equations'><span class="jlbinding">HarmonicBalance.get_harmonic_equations</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_harmonic_equations(diff_eom::DifferentialEquation; fast_time=nothing, slow_time=nothing)
```


Apply the harmonic ansatz, followed by the slow-flow, Fourier transform and dropping higher-order derivatives to obtain a set of ODEs (the harmonic equations) governing the harmonics of `diff_eom`.

The harmonics evolve in `slow_time`, the oscillating terms themselves in `fast_time`. If no input is used, a variable T is defined for `slow_time` and `fast_time` is taken as the independent variable of `diff_eom`.

By default, all products of order &gt; 1 of `slow_time`-derivatives are dropped, which means the equations are linear in the time-derivatives.

**Example**

```julia
julia> @variables t, x(t), ω0, ω, F;

# enter the simple harmonic oscillator
julia> diff_eom = DifferentialEquation( d(x,t,2) + ω0^2 * x ~ F *cos(ω*t), x);

# expand x in the harmonic ω
julia> add_harmonic!(diff_eom, x, ω);

# get equations for the harmonics evolving in the slow time T
julia> harmonic_eom = get_harmonic_equations(diff_eom)

A set of 2 harmonic equations
Variables: u1(T), v1(T)
Parameters: ω0, ω, F

Harmonic ansatz:
x(t) = u1*cos(ωt) + v1*sin(ωt)

Harmonic equations:

(ω0^2)*u1(T) + (2//1)*ω*Differential(T)(v1(T)) - (ω^2)*u1(T) ~ F

(ω0^2)*v1(T) - (ω^2)*v1(T) - (2//1)*ω*Differential(T)(u1(T)) ~ 0
```



[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicEquation.jl#L269-L310)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.KrylovBogoliubov.get_krylov_equations' href='#HarmonicBalance.KrylovBogoliubov.get_krylov_equations'><span class="jlbinding">HarmonicBalance.KrylovBogoliubov.get_krylov_equations</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_krylov_equations(
    diff_eom::DifferentialEquation;
    order,
    fast_time,
    slow_time
)

```


Apply the Krylov-Bogoliubov averaging method to a specific `order` to obtain a set of ODEs (the slow-flow equations) governing the harmonics of `diff_eom`.

The harmonics evolve in `slow_time`, the oscillating terms themselves in `fast_time`. If no input is used, a variable T is defined for `slow_time` and `fast_time` is taken as the independent variable of `diff_eom`.

Krylov-Bogoliubov averaging method can be applied up to `order = 2`.

**Example**

```julia
julia> @variables t, x(t), ω0, ω, F;

# enter the simple harmonic oscillator
julia> diff_eom = DifferentialEquation( d(x,t,2) + ω0^2 * x ~ F *cos(ω*t), x);

# expand x in the harmonic ω
julia> add_harmonic!(diff_eom, x, ω);

# get equations for the harmonics evolving in the slow time T to first order
julia> harmonic_eom = get_krylov_equations(diff_eom, order = 1)

A set of 2 harmonic equations
Variables: u1(T), v1(T)
Parameters: ω, F, ω0

Harmonic ansatz:
xˍt(t) =
x(t) = u1(T)*cos(ωt) + v1(T)*sin(ωt)

Harmonic equations:

((1//2)*(ω^2)*v1(T) - (1//2)*(ω0^2)*v1(T)) / ω ~ Differential(T)(u1(T))

((1//2)*(ω0^2)*u1(T) - (1//2)*F - (1//2)*(ω^2)*u1(T)) / ω ~ Differential(T)(v1(T))
```



[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/KrylovBogoliubov.jl#L25)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.add_harmonic!' href='#HarmonicBalance.add_harmonic!'><span class="jlbinding">HarmonicBalance.add_harmonic!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
add_harmonic!(diff_eom::DifferentialEquation, var::Num, ω)

```


Add the harmonic `ω` to the harmonic ansatz used to expand the variable `var` in `diff_eom`.

**Example**

**define the simple harmonic oscillator and specify that x(t) oscillates with frequency ω**

```julia
julia> @variables t, x(t), y(t), ω0, ω, F, k;
julia> diff_eq = DifferentialEquation(d(x,t,2) + ω0^2 * x ~ F * cos(ω*t), x);
julia> add_harmonic!(diff_eq, x, ��) # expand x using ω

System of 1 differential equations
Variables:       x(t)
Harmonic ansatz: x(t) => ω;

(ω0^2)*x(t) + Differential(t)(Differential(t)(x(t))) ~ F*cos(t*ω)
```



[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/DifferentialEquation.jl#L88)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.get_independent_variables' href='#HarmonicBalance.get_independent_variables'><span class="jlbinding">HarmonicBalance.get_independent_variables</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_independent_variables(
    diff_eom::DifferentialEquation
) -> Any

```


Return the independent dependent variables of `diff_eom`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/DifferentialEquation.jl#L131)



```julia
get_independent_variables(
    eom::HarmonicEquation
) -> Vector{Num}

```


Return the independent variables (typically time) of `eom`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicEquation.jl#L221)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Symbolics.get_variables' href='#Symbolics.get_variables'><span class="jlbinding">Symbolics.get_variables</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_variables(diff_eom::DifferentialEquation) -> Vector{Num}

```


Return the dependent variables of `diff_eom`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/DifferentialEquation.jl#L113)



```julia
get_variables(eom::HarmonicEquation) -> Vector{Num}

```


Get the internal symbols of the independent variables of `eom`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicEquation.jl#L188)

</details>


## Solving and transforming solutions {#Solving-and-transforming-solutions}
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.get_steady_states' href='#HarmonicBalance.get_steady_states'><span class="jlbinding">HarmonicBalance.get_steady_states</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_steady_states(problem::HarmonicEquation,
                    method::HarmonicBalanceMethod,
                    swept_parameters,
                    fixed_parameters;
                    show_progress=true,
                    sorting="nearest",
                    classify_default=true)
```


Solves `problem` with the `method` over the ranges specified by `swept_parameters`, keeping `fixed_parameters` constant. `swept_parameters` accepts pairs mapping symbolic variables to arrays or ranges. `fixed_parameters` accepts pairs mapping symbolic variables to numbers.

**Keyword arguments**
- `show_progress`: Indicate whether a progress bar should be displayed.
  
- `sorting`: the method used by `sort_solutions` to get continuous solutions branches.   The current options are `"hilbert"` (1D sorting along a Hilbert curve), `"nearest"`   (nearest-neighbor sorting) and `"none"`.
  
- `classify_default`: If `true`, the solutions will be classified using the default   classification method.
  

**Example**

solving a simple harmonic oscillator $m \ddot{x} + γ \dot{x} + ω_0^2 x = F \cos(ωt)$ to obtain the response as a function of $ω$

```julia
# having obtained a Problem object, let's find steady states
julia> range = (ω => range(0.8, 1.2, 100) ) # 100 parameter sets to solve
julia> fixed = ParameterList(m => 1, γ => 0.01, F => 0.5, ω_0 => 1)
julia> get_steady_states(problem, range, fixed)

A steady state result for 100 parameter points

    Solution branches:   1
       of which real:    1
       of which stable:  1

    Classes: stable, physical, Hopf, binary_labels

```


It is also possible to perform 2-dimensional sweeps.

```julia
# The swept parameters take precedence over fixed -> use the same fixed
julia> range = (ω => range(0.8,1.2,100), F => range(0.1,1.0,10) )

# The swept parameters take precedence over fixed -> the F in fixed is now ignored
julia> get_steady_states(problem, range, fixed)

A steady state result for 1000 parameter points

    Solution branches:   1
       of which real:    1
       of which stable:  1

    Classes: stable, physical, Hopf, binary_labels
```



[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/solve_homotopy.jl#L1-L60)

</details>


### Methods
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.WarmUp' href='#HarmonicBalance.WarmUp'><span class="jlbinding">HarmonicBalance.WarmUp</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.TotalDegree' href='#HarmonicBalance.TotalDegree'><span class="jlbinding">HarmonicBalance.TotalDegree</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.Polyhedral' href='#HarmonicBalance.Polyhedral'><span class="jlbinding">HarmonicBalance.Polyhedral</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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


### Access solutions {#Access-solutions}
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.get_single_solution' href='#HarmonicBalance.get_single_solution'><span class="jlbinding">HarmonicBalance.get_single_solution</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_single_solution(
    res::HarmonicBalance.Result{S, P, D, F} where {D, F<:FunctionWrappers.FunctionWrapper{Array{S, 2}, Tuple{Array{S, 1}}}};
    branch,
    index
)

```


Return an ordered dictionary specifying all variables and parameters of the solution in `result` on `branch` at the position `index`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/transform_solutions.jl#L1)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.transform_solutions' href='#HarmonicBalance.transform_solutions'><span class="jlbinding">HarmonicBalance.transform_solutions</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
transform_solutions(
    res::HarmonicBalance.Result{S, ParType, D, F} where {ParType<:Number, D, F<:FunctionWrappers.FunctionWrapper{Array{S, 2}, Tuple{Array{S, 1}}}},
    func;
    branches,
    realify
) -> Vector

```


Takes a `Result` object and a string `f` representing a Symbolics.jl expression. Returns an array with the values of `f` evaluated for the respective solutions. Additional substitution rules can be specified in `rules` in the format `("a" => val)` or `(a => val)`


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/transform_solutions.jl#L64)

</details>


### Classify
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.classify_solutions!' href='#HarmonicBalance.classify_solutions!'><span class="jlbinding">HarmonicBalance.classify_solutions!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
classify_solutions!(
    res::HarmonicBalance.Result,
    func::Union{Function, String},
    name::String;
    physical
) -> Any

```


Creates a solution class in `res` using the function `func` (parsed into Symbolics.jl input). The new class is labeled with `name` and stored under `res.classes[name]`. By default, only physical (real) solutions are classified, and `false` is returned for the rest. To also classify complex solutions, set `physical=false`.

**Example**

```julia
# solve a previously-defined problem
res = get_steady_states(problem, swept_parameters, fixed_parameters)

# classify, store in result.classes["large_amplitude"]
classify_solutions!(res, "sqrt(u1^2 + v1^2) > 1.0" , "large_amplitude")
```



[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/classification.jl#L1)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.get_class' href='#HarmonicBalance.get_class'><span class="jlbinding">HarmonicBalance.get_class</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_class(
    res::HarmonicBalance.Result,
    branch::Int64,
    class::String
) -> Any

```


Returns an array of booleans classifying `branch` in the solutions in `res` according to `class`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/classification.jl#L43)



```julia
get_class(
    soln::HarmonicBalance.Result,
    class::String
) -> Vector

```


Returns an array of booleans classifying each branch in the solutions in `res` according to `class`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/classification.jl#L53)

</details>


### Plotting
<details class='jldocstring custom-block' open>
<summary><a id='RecipesBase.plot' href='#RecipesBase.plot'><span class="jlbinding">RecipesBase.plot</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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



```julia
plot(
    res::HarmonicBalance.Result{S, P, D, F} where F<:FunctionWrappers.FunctionWrapper{Array{S, 2}, Tuple{Array{S, 1}}},
    varargs...;
    cut,
    kwargs...
) -> Plots.Plot

```


**Plot a `Result` object.**

Class selection done by passing `String` or `Vector{String}` as kwarg:

```
class       :   only plot solutions in this class(es) ("all" --> plot everything)
not_class   :   do not plot solutions in this class(es)
```


Other kwargs are passed onto Plots.gr().

See also `plot!`

The x,y,z arguments are Strings compatible with Symbolics.jl, e.g., `y=2*sqrt(u1^2+v1^2)` plots the amplitude of the first quadratures multiplied by 2.

**1D plots**

```
plot(res::Result; x::String, y::String, class="default", not_class=[], kwargs...)
plot(res::Result, y::String; kwargs...) # take x automatically from Result
```


Default behaviour is to plot stable solutions as full lines, unstable as dashed.

If a sweep in two parameters were done, i.e., `dim(res)==2`, a one dimensional cut can be plotted by using the keyword `cut` were it takes a `Pair{Num, Float}` type entry. For example, `plot(res, y="sqrt(u1^2+v1^2), cut=(λ => 0.2))` plots a cut at `λ = 0.2`.

****

**2D plots**

```
plot(res::Result; z::String, branch::Int64, class="physical", not_class=[], kwargs...)
```


To make the 2d plot less chaotic it is required to specify the specific `branch` to plot, labeled by a `Int64`.

The x and y axes are taken automatically from `res`


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/plotting_Plots.jl#L11)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RecipesBase.plot!' href='#RecipesBase.plot!'><span class="jlbinding">RecipesBase.plot!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
plot!(
    res::HarmonicBalance.Result,
    varargs...;
    kwargs...
) -> Plots.Plot

```


Similar to `plot` but adds a plot onto an existing plot.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/plotting_Plots.jl#L62)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.plot_phase_diagram' href='#HarmonicBalance.plot_phase_diagram'><span class="jlbinding">HarmonicBalance.plot_phase_diagram</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
plot_phase_diagram(
    res::HarmonicBalance.Result;
    kwargs...
) -> Plots.Plot

```


Plot the number of solutions in a `Result` object as a function of the parameters. Works with 1D and 2D datasets.

Class selection done by passing `String` or `Vector{String}` as kwarg:

```
class::String       :   only count solutions in this class ("all" --> plot everything)
not_class::String   :   do not count solutions in this class
```


Other kwargs are passed onto Plots.gr()


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/plotting_Plots.jl#L275)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.plot_spaghetti' href='#HarmonicBalance.plot_spaghetti'><span class="jlbinding">HarmonicBalance.plot_spaghetti</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
plot_spaghetti(res::Result; x, y, z, kwargs...)
```


Plot a three dimension line plot of a `Result` object as a function of the parameters. Works with 1D and 2D datasets.

Class selection done by passing `String` or `Vector{String}` as kwarg:

```
class::String       :   only count solutions in this class ("all" --> plot everything)
not_class::String   :   do not count solutions in this class
```


Other kwargs are passed onto Plots.gr()


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/plotting_Plots.jl#L342-L354)

</details>


## Limit cycles {#Limit-cycles}
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LimitCycles.get_limit_cycles' href='#HarmonicBalance.LimitCycles.get_limit_cycles'><span class="jlbinding">HarmonicBalance.LimitCycles.get_limit_cycles</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_limit_cycles(
    eom::HarmonicEquation, method::HarmonicBalanceMethod, swept, fixed, ω_lc; kwargs...)
```


Variant of `get_steady_states` for a limit cycle problem characterised by a Hopf frequency (usually called ω_lc)

Solutions with ω_lc = 0 are labelled unphysical since this contradicts the assumption of distinct harmonic variables corresponding to distinct harmonics.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LimitCycles/gauge_fixing.jl#L92-L99)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LimitCycles.get_cycle_variables' href='#HarmonicBalance.LimitCycles.get_cycle_variables'><span class="jlbinding">HarmonicBalance.LimitCycles.get_cycle_variables</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_cycle_variables(
    eom::HarmonicEquation,
    ω_lc::Num
) -> Vector{HarmonicVariable}

```


Return the harmonic variables which participate in the limit cycle labelled by `ω_lc`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LimitCycles/gauge_fixing.jl#L22)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LimitCycles.add_pairs!' href='#HarmonicBalance.LimitCycles.add_pairs!'><span class="jlbinding">HarmonicBalance.LimitCycles.add_pairs!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
add_pairs!(eom::DifferentialEquation; ω_lc::Num, n=1)
```


Add a limit cycle harmonic `ω_lc` to the system Equivalent to adding `n` pairs of harmonics ω +- ω_lc for each existing ω.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LimitCycles/gauge_fixing.jl#L9-L14)

</details>


## Linear Response {#Linear-Response}
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LinearResponse.plot_eigenvalues-Union{Tuple{HarmonicBalance.Result{S, P, D, F} where {D, F<:FunctionWrappers.FunctionWrapper{Matrix{S}, Tuple{Vector{S}}}}}, Tuple{P}, Tuple{S}} where {S, P}' href='#HarmonicBalance.LinearResponse.plot_eigenvalues-Union{Tuple{HarmonicBalance.Result{S, P, D, F} where {D, F<:FunctionWrappers.FunctionWrapper{Matrix{S}, Tuple{Vector{S}}}}}, Tuple{P}, Tuple{S}} where {S, P}'><span class="jlbinding">HarmonicBalance.LinearResponse.plot_eigenvalues</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_eigenvalues(res::Result; branch::Int, class=["physical"], type=:imag, projection=v -> 1, cscheme=:default, kwargs...)
```


Plot the eigenvalues of the jacobian in the rotating frame for Result `res` on `branch`. Either the real (`type=:real``) or imaginary part (`type=:imag``) can be plotted. The`projection` function ℜᵈ → ℜ is applied to the eigenvectors and defines the color of the eigenvalues. The color scheme can be set to a custom one or to the default one.

Any kwargs are fed to Plots&#39; gr().

Solutions not belonging to the `physical` class are ignored.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LinearResponse/plotting.jl#L259-L267)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LinearResponse.plot_linear_response-Tuple{HarmonicBalance.Result, Num}' href='#HarmonicBalance.LinearResponse.plot_linear_response-Tuple{HarmonicBalance.Result, Num}'><span class="jlbinding">HarmonicBalance.LinearResponse.plot_linear_response</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_linear_response(res::Result, nat_var::Num; Ω_range, branch::Int, order=1, logscale=false, show_progress=true, kwargs...)
```


Plot the linear response to white noise of the variable `nat_var` for Result `res` on `branch` for input frequencies `Ω_range`. Slow-time derivatives up to `order` are kept in the process.

Any kwargs are fed to Plots&#39; gr().

Solutions not belonging to the `physical` class are ignored.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LinearResponse/plotting.jl#L120-L129)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LinearResponse.plot_rotframe_jacobian_response-Union{Tuple{HarmonicBalance.Result{S, P, D, F} where {D, F<:FunctionWrappers.FunctionWrapper{Matrix{S}, Tuple{Vector{S}}}}}, Tuple{P}, Tuple{S}} where {S, P}' href='#HarmonicBalance.LinearResponse.plot_rotframe_jacobian_response-Union{Tuple{HarmonicBalance.Result{S, P, D, F} where {D, F<:FunctionWrappers.FunctionWrapper{Matrix{S}, Tuple{Vector{S}}}}}, Tuple{P}, Tuple{S}} where {S, P}'><span class="jlbinding">HarmonicBalance.LinearResponse.plot_rotframe_jacobian_response</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_rotframe_jacobian_response(res::Result; Ω_range, branch::Int, logscale=true, damping_mod = 1.0, show_progress=true, kwargs...)
```


Plot the linear response to white noise in the rotating frame for Result `res` on `branch` for input frequencies `Ω_range`. &#39;damping_mod&#39; gets multiplied by the real part of the eigenvalues of the Jacobian in order to be able to make peaks with similar frequency separately identifiable.

Any kwargs are fed to Plots&#39; gr().

Solutions not belonging to the `physical` class are ignored.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LinearResponse/plotting.jl#L214-L222)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.get_Jacobian' href='#HarmonicBalance.get_Jacobian'><span class="jlbinding">HarmonicBalance.get_Jacobian</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_Jacobian(eom)

```


Obtain the symbolic Jacobian matrix of `eom`. This is the linearised left-hand side of F(u) = du/dT.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/Jacobian.jl#L46)



Obtain a Jacobian from a `DifferentialEquation` by first converting it into a `HarmonicEquation`. 


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/Jacobian.jl#L60)



Get the Jacobian of a set of equations `eqs` with respect to the variables `vars`. 


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/Jacobian.jl#L69)

</details>


## Extensions

### OrdinaryDiffEq
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.AdiabaticSweep' href='#HarmonicBalance.AdiabaticSweep'><span class="jlbinding">HarmonicBalance.AdiabaticSweep</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.follow_branch' href='#HarmonicBalance.follow_branch'><span class="jlbinding">HarmonicBalance.follow_branch</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
follow_branch(
    starting_branch::Int64,
    res::HarmonicBalance.Result;
    y,
    sweep,
    tf,
    ϵ
) -> Tuple{Any, Any}

```


Return the indexes and values following stable branches along a 1D sweep. When a no stable solutions are found (e.g. in a bifurcation), the next stable solution is calculated by time evolving the previous solution (quench).

**Keyword arguments**
- `y`:  Dependent variable expression (parsed into Symbolics.jl) to evaluate the followed solution branches on .
  
- `sweep`: Direction for the sweeping of solutions. A `right` (`left`) sweep proceeds from the first (last) solution, ordered as the sweeping parameter.
  
- `tf`: time to reach steady
  
- `ϵ`: small random perturbation applied to quenched solution, in a bifurcation in order to favour convergence in cases where multiple solutions are identically accessible (e.g. symmetry breaking into two equal amplitude states)
  


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/TimeEvolution/hysteresis_sweep.jl#L15)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.plot_1D_solutions_branch' href='#HarmonicBalance.plot_1D_solutions_branch'><span class="jlbinding">HarmonicBalance.plot_1D_solutions_branch</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



1D plot with the followed branch highlighted


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/TimeEvolution/plotting.jl#L62-L64)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.follow_branch-Tuple{Int64, HarmonicBalance.Result}-manual-API' href='#HarmonicBalance.follow_branch-Tuple{Int64, HarmonicBalance.Result}-manual-API'><span class="jlbinding">HarmonicBalance.follow_branch</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
follow_branch(
    starting_branch::Int64,
    res::HarmonicBalance.Result;
    y,
    sweep,
    tf,
    ϵ
) -> Tuple{Any, Any}

```


Return the indexes and values following stable branches along a 1D sweep. When a no stable solutions are found (e.g. in a bifurcation), the next stable solution is calculated by time evolving the previous solution (quench).

**Keyword arguments**
- `y`:  Dependent variable expression (parsed into Symbolics.jl) to evaluate the followed solution branches on .
  
- `sweep`: Direction for the sweeping of solutions. A `right` (`left`) sweep proceeds from the first (last) solution, ordered as the sweeping parameter.
  
- `tf`: time to reach steady
  
- `ϵ`: small random perturbation applied to quenched solution, in a bifurcation in order to favour convergence in cases where multiple solutions are identically accessible (e.g. symmetry breaking into two equal amplitude states)
  


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/TimeEvolution/hysteresis_sweep.jl#L15)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='RecipesBase.plot-Tuple{ODESolution, Any, HarmonicEquation}-manual-API' href='#RecipesBase.plot-Tuple{ODESolution, Any, HarmonicEquation}-manual-API'><span class="jlbinding">RecipesBase.plot</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



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


### SteadyStateSweep
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.steady_state_sweep' href='#HarmonicBalance.steady_state_sweep'><span class="jlbinding">HarmonicBalance.steady_state_sweep</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
steady_state_sweep(prob::SteadyStateProblem, alg::DynamicSS; varied::Pair, kwargs...)
```


Sweeps through a range of parameter values using a dynamic steady state solver `DynamicSS` of the `SteadyStateDiffEq.jl` package. Given a steady state problem and a parameter to vary, computes the steady state solution for each value in the sweep range. The solutions are returned as a vector where each element corresponds to the steady state found at that parameter value.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/SteadyStateDiffEqExt.jl#L10-L18)



```julia
steady_state_sweep(prob_np::NonlinearProblem, prob_ss::SteadyStateProblem,
                  alg_np, alg_ss::DynamicSS; varied::Pair, kwargs...)
```


Performs a parameter sweep by combining nonlinear root `alg_np` and steady state solvers `alg_ss`. For each parameter value, it first attempts a direct nonlinear root solver and checks its stability. If the solution is unstable or not found, it switches to a dynamic steady state solver. This hybrid approach is much faster then only using a steady state solver.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/SteadyStateDiffEqExt.jl#L37-L45)

</details>


### ModelingToolkit
<details class='jldocstring custom-block' open>
<summary><a id='SciMLBase.ODEProblem' href='#SciMLBase.ODEProblem'><span class="jlbinding">SciMLBase.ODEProblem</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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



```julia
ODEProblem(
    eom::Union{DifferentialEquation, HarmonicEquation},
    u0,
    tspan::Tuple,
    p::AbstractDict;
    in_place,
    kwargs...
) -> Any

```


Creates and ModelingToolkit.ODEProblem from a DifferentialEquation.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/ModelingToolkitExt.jl#L100)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ModelingToolkit.ODESystem' href='#ModelingToolkit.ODESystem'><span class="jlbinding">ModelingToolkit.ODESystem</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ODESystem(eom::HarmonicEquation) -> Any

```


Creates and ModelingToolkit.ODESystem from a HarmonicEquation.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/ModelingToolkitExt.jl#L39)



```julia
ODESystem(diff_eq::DifferentialEquation) -> Any

```


Creates and ModelingToolkit.ODESystem from a DifferentialEquation.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/ModelingToolkitExt.jl#L69)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SciMLBase.SteadyStateProblem' href='#SciMLBase.SteadyStateProblem'><span class="jlbinding">SciMLBase.SteadyStateProblem</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SteadyStateProblem(
    eom::HarmonicEquation,
    u0,
    p::AbstractDict;
    in_place,
    kwargs...
) -> Any

```


Creates and ModelingToolkit.SteadyStateProblem from a DifferentialEquation.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/ModelingToolkitExt.jl#L135)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='SciMLBase.NonlinearProblem' href='#SciMLBase.NonlinearProblem'><span class="jlbinding">SciMLBase.NonlinearProblem</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
NonlinearProblem(
    eom::HarmonicEquation,
    u0,
    p::AbstractDict;
    in_place,
    kwargs...
) -> Any

```


Creates and ModelingToolkit.NonlinearProblem from a DifferentialEquation.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/ext/ModelingToolkitExt.jl#L123)

</details>

