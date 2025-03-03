
# Linear response (WIP) {#linresp_man}

This module currently has two goals. One is calculating the first-order Jacobian, used to obtain stability and approximate (but inexpensive) the linear response of steady states. The other is calculating the full response matrix as a function of frequency; this is more accurate but more expensive. 

The methodology used is explained in [Jan Kosata phd thesis](https://www.doi.org/10.3929/ethz-b-000589190).

## Stability

The Jacobian is used to evaluate stability of the solutions. It can be shown explicitly,
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.get_Jacobian-manual-linear_response' href='#HarmonicBalance.get_Jacobian-manual-linear_response'><span class="jlbinding">HarmonicBalance.get_Jacobian</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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


## Linear response {#Linear-response}

The response to white noise can be shown with `plot_linear_response`. Depending on the `order` argument, different methods are used. 
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LinearResponse.plot_linear_response-manual-linear_response' href='#HarmonicBalance.LinearResponse.plot_linear_response-manual-linear_response'><span class="jlbinding">HarmonicBalance.LinearResponse.plot_linear_response</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
plot_linear_response(res::Result, nat_var::Num; Ω_range, branch::Int, order=1, logscale=false, show_progress=true, kwargs...)
```


Plot the linear response to white noise of the variable `nat_var` for Result `res` on `branch` for input frequencies `Ω_range`. Slow-time derivatives up to `order` are kept in the process.

Any kwargs are fed to Plots&#39; gr().

Solutions not belonging to the `physical` class are ignored.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LinearResponse/plotting.jl#L120-L129)

</details>


### First order {#First-order}

The simplest way to extract the linear response of a steady state is to evaluate the Jacobian of the harmonic equations. Each of its eigenvalues $\lambda$ describes a Lorentzian peak in the response; $\text{Re}[\lambda]$ gives its center and $\text{Im}[\lambda]$ its width. Transforming the harmonic variables into the non-rotating frame (that is, inverting the harmonic ansatz) then gives the response as it would be observed in an experiment.

The advantage of this method is that for a given parameter set, only one matrix diagonalization is needed to fully describe the response spectrum. However, the method is inaccurate for response frequencies far from the frequencies used in the harmonic ansatz (it relies on the response oscillating slowly in the rotating frame). 

Behind the scenes, the spectra are stored using the dedicated structs `Lorentzian` and `JacobianSpectrum`.
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LinearResponse.JacobianSpectrum-manual-linear_response' href='#HarmonicBalance.LinearResponse.JacobianSpectrum-manual-linear_response'><span class="jlbinding">HarmonicBalance.LinearResponse.JacobianSpectrum</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
mutable struct JacobianSpectrum{T<:Real}
```


Holds a set of `Lorentzian` objects belonging to a variable.

**Fields**
- `peaks::Array{HarmonicBalance.LinearResponse.Lorentzian{T}, 1} where T<:Real`
  

**Constructor**

```julia
JacobianSpectrum(res::Result; index::Int, branch::Int)
```



[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LinearResponse/types.jl#L21)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LinearResponse.Lorentzian-manual-linear_response' href='#HarmonicBalance.LinearResponse.Lorentzian-manual-linear_response'><span class="jlbinding">HarmonicBalance.LinearResponse.Lorentzian</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct Lorentzian{T<:Real}
```


Holds the three parameters of a Lorentzian peak, defined as A / sqrt((ω-ω0)² + Γ²).

**Fields**
- `ω0::Real`
  
- `Γ::Real`
  
- `A::Real`
  


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LinearResponse/types.jl#L1)

</details>


### Higher orders {#Higher-orders}

Setting `order > 1` increases the accuracy of the response spectra. However, unlike for the Jacobian, here we must perform a matrix inversion for each response frequency.  
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LinearResponse.ResponseMatrix-manual-linear_response' href='#HarmonicBalance.LinearResponse.ResponseMatrix-manual-linear_response'><span class="jlbinding">HarmonicBalance.LinearResponse.ResponseMatrix</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct ResponseMatrix
```


Holds the compiled response matrix of a system.

**Fields**
- `matrix::Matrix{Function}`: The response matrix (compiled).
  
- `symbols::Vector{Num}`: Any symbolic variables in `matrix` to be substituted at evaluation.
  
- `variables::Vector{HarmonicVariable}`: The frequencies of the harmonic variables underlying `matrix`. These are needed to transform the harmonic variables to the non-rotating frame.
  


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LinearResponse/types.jl#L40)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LinearResponse.get_response-manual-linear_response' href='#HarmonicBalance.LinearResponse.get_response-manual-linear_response'><span class="jlbinding">HarmonicBalance.LinearResponse.get_response</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_response(
    rmat::HarmonicBalance.LinearResponse.ResponseMatrix,
    s::OrderedCollections.OrderedDict,
    Ω
) -> Any

```


For `rmat` and a solution dictionary `s`, calculate the total response to a perturbative force at frequency `Ω`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LinearResponse/response.jl#L63)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.LinearResponse.get_response_matrix-manual-linear_response' href='#HarmonicBalance.LinearResponse.get_response_matrix-manual-linear_response'><span class="jlbinding">HarmonicBalance.LinearResponse.get_response_matrix</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_response_matrix(diff_eq::DifferentialEquation, freq::Num; order=2)
```


Obtain the symbolic linear response matrix of a `diff_eq` corresponding to a perturbation frequency `freq`. This routine cannot accept a `HarmonicEquation` since there, some time-derivatives are already dropped. `order` denotes the highest differential order to be considered.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/LinearResponse/response.jl#L1-L8)

</details>

