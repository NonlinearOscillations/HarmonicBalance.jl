
# Entering equations of motion {#Entering-equations-of-motion}

The struct `DifferentialEquation` is the primary input method; it holds an ODE or a coupled system of ODEs composed of terms with harmonic time-dependence The dependent variables are specified during input, any other symbols are identified as parameters. Information on which variable is to be expanded in which harmonic is specified using `add_harmonic!`.

`DifferentialEquation.equations` stores a dictionary assigning variables to equations. This information is necessary because the harmonics belonging to a variable are later used to Fourier-transform its corresponding ODE.
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.d-manual-entering_eom' href='#HarmonicBalance.d-manual-entering_eom'><span class="jlbinding">HarmonicBalance.d</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



The derivative of f w.r.t. x of degree deg


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/HarmonicVariable.jl#L108)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.DifferentialEquation-manual-entering_eom' href='#HarmonicBalance.DifferentialEquation-manual-entering_eom'><span class="jlbinding">HarmonicBalance.DifferentialEquation</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



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
<summary><a id='HarmonicBalance.add_harmonic!-manual-entering_eom' href='#HarmonicBalance.add_harmonic!-manual-entering_eom'><span class="jlbinding">HarmonicBalance.add_harmonic!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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
<summary><a id='Symbolics.get_variables-Tuple{DifferentialEquation}-manual-entering_eom' href='#Symbolics.get_variables-Tuple{DifferentialEquation}-manual-entering_eom'><span class="jlbinding">Symbolics.get_variables</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_variables(diff_eom::DifferentialEquation) -> Vector{Num}

```


Return the dependent variables of `diff_eom`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/DifferentialEquation.jl#L113)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.get_independent_variables-Tuple{DifferentialEquation}-manual-entering_eom' href='#HarmonicBalance.get_independent_variables-Tuple{DifferentialEquation}-manual-entering_eom'><span class="jlbinding">HarmonicBalance.get_independent_variables</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_independent_variables(
    diff_eom::DifferentialEquation
) -> Any

```


Return the independent dependent variables of `diff_eom`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/DifferentialEquation.jl#L131)

</details>

