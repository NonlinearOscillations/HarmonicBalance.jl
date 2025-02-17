
# Extracting harmonic equations {#Extracting-harmonic-equations}

## Harmonic Balance method {#Harmonic-Balance-method}

Once a `DifferentialEquation` is defined and its harmonics specified, one can extract the harmonic equations using `get_harmonic_equations`, which itself is composed of the subroutines `harmonic_ansatz`, `slow_flow`, `fourier_transform!` and `drop_powers`. 

The harmonic equations use an additional time variable specified as `slow_time` in `get_harmonic_equations`. This is essentially a label distinguishing the time dependence of the harmonic variables (expected to be slow) from that of the oscillating terms (expected to be fast). When the equations are Fourier-transformed to remove oscillating terms, `slow_time` is treated as a constant. Such an approach is exact when looking for steady states. 
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.get_harmonic_equations-manual-extracting_harmonics' href='#HarmonicBalance.get_harmonic_equations-manual-extracting_harmonics'><span class="jlbinding">HarmonicBalance.get_harmonic_equations</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
<summary><a id='HarmonicBalance.harmonic_ansatz-manual-extracting_harmonics' href='#HarmonicBalance.harmonic_ansatz-manual-extracting_harmonics'><span class="jlbinding">HarmonicBalance.harmonic_ansatz</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
harmonic_ansatz(eom::DifferentialEquation, time::Num; coordinates="Cartesian")
```


Expand each variable of `diff_eom` using the harmonics assigned to it with `time` as the time variable. For each harmonic of each variable, instance(s) of `HarmonicVariable` are automatically created and named.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicEquation.jl#L56-L62)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.slow_flow-manual-extracting_harmonics' href='#HarmonicBalance.slow_flow-manual-extracting_harmonics'><span class="jlbinding">HarmonicBalance.slow_flow</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
slow_flow(eom::HarmonicEquation; fast_time::Num, slow_time::Num, degree=2)
```


Removes all derivatives w.r.t `fast_time` (and their products) in `eom` of power `degree`. In the remaining derivatives, `fast_time` is replaced by `slow_time`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicEquation.jl#L123-L128)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.fourier_transform-manual-extracting_harmonics' href='#HarmonicBalance.fourier_transform-manual-extracting_harmonics'><span class="jlbinding">HarmonicBalance.fourier_transform</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fourier_transform(
    eom::HarmonicEquation,
    time::Num
) -> HarmonicEquation

```


Extract the Fourier components of `eom` corresponding to the harmonics specified in `eom.variables`. For each non-zero harmonic of each variable, 2 equations are generated (cos and sin Fourier coefficients). For each zero (constant) harmonic, 1 equation is generated `time` does not appear in the resulting equations anymore.

Underlying assumption: all time-dependences are harmonic.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicEquation.jl#L230)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.ExprUtils.drop_powers-manual-extracting_harmonics' href='#HarmonicBalance.ExprUtils.drop_powers-manual-extracting_harmonics'><span class="jlbinding">HarmonicBalance.ExprUtils.drop_powers</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
drop_powers(expr, vars, deg)

```


Remove parts of `expr` where the combined power of `vars` is =&gt; `deg`.

**Example**

```julia
julia> @variables x,y;
julia>drop_powers((x+y)^2, x, 2)
y^2 + 2*x*y
julia>drop_powers((x+y)^2, [x,y], 2)
0
julia>drop_powers((x+y)^2 + (x+y)^3, [x,y], 3)
x^2 + y^2 + 2*x*y
```



[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/ExprUtils/drop_powers.jl#L1)

</details>


## HarmonicVariable and HarmonicEquation types {#HarmonicVariable-and-HarmonicEquation-types}

The equations governing the harmonics are stored using the two following structs. When going from the original to the harmonic equations, the harmonic ansatz $x_i(t) = \sum_{j=1}^M u_{i,j}  (T)  \cos(\omega_{i,j} t)+ v_{i,j}(T) \sin(\omega_{i,j} t)$ is used. Internally, each pair $(u_{i,j}, v_{i,j})$ is stored as a `HarmonicVariable`. This includes the identification of $\omega_{i,j}$ and $x_i(t)$, which is needed to later reconstruct $x_i(t)$.
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.HarmonicVariable-manual-extracting_harmonics' href='#HarmonicBalance.HarmonicVariable-manual-extracting_harmonics'><span class="jlbinding">HarmonicBalance.HarmonicVariable</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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


When the full set of equations of motion is expanded using the harmonic ansatz, the result is stored as a `HarmonicEquation`. For an initial equation of motion consisting of $M$ variables, each expanded in $N$ harmonics, the resulting `HarmonicEquation` holds $2NM$ equations of $2NM$ variables. Each symbol not corresponding to a variable is identified as a parameter.

A `HarmonicEquation` can be either parsed into a steady-state `Problem` or solved using a dynamical ODE solver.
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.HarmonicEquation-manual-extracting_harmonics' href='#HarmonicBalance.HarmonicEquation-manual-extracting_harmonics'><span class="jlbinding">HarmonicBalance.HarmonicEquation</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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



MarkdownAST.LineBreak()


