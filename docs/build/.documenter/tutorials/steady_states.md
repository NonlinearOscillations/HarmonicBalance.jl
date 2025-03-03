
# Finding the staedy states of a Duffing oscillator {#Duffing}

Here we show the workflow of HarmonicBalance.jl on a simple example - the driven Duffing oscillator. The equation of motion for the displacement $x(t)$ reads

$$\begin{equation}
\underbrace{\ddot{x}(t) + \gamma \dot{x}(t) + \omega_0^2 x(t)}_{\text{damped harmonic oscillator}} + \underbrace{\alpha x(t)^3}_{\text{Duffing coefficient}} = \underbrace{F \cos(\omega t)}_{\text{periodic drive}}
\end{equation}$$

In general, there is no analytical solution to the differential equation. Fortunately, some harmonics are more important than others. By truncating the infinite-dimensional Fourier space to a set of judiciously chosen harmonics, we may obtain a soluble system. For the Duffing resonator, we can well try to only consider the drive frequency $\omega$. To implement this, we use the _harmonic ansatz_

$$x(t) = U \cos(\omega t) + V \sin(\omega t) \,,$$

which constraints the spectrum of $x(t)$ to a single harmonic. Fixing the quadratures $U$ and $V$ to be constant then reduces the differential equation to two coupled cubic polynomial equations (for more details on this step, see the appendices in the [white paper](https://scipost.org/SciPostPhysCodeb.6)). Finding the roots of coupled polynomials is in general very hard. We here apply the method of homotopy continuation, as implemented in [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/) which is guaranteed to find the complete set of roots.

First we need to declare the symbolic variables (the excellent [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) is used here).

```julia
using HarmonicBalance
@variables α ω ω0 F γ t x(t) # declare constant variables and a function x(t)
```


```
<< @example-block not executed in draft mode >>
```


Next, we have to input the equations of motion. This will be stored as a `DifferentialEquation`. The input needs to specify that only `x` is a mathematical variable, the other symbols are parameters:

```julia
diff_eq = DifferentialEquation(d(x,t,2) + ω0^2*x + α*x^3 + γ*d(x,t) ~ F*cos(ω*t), x)
```


```
<< @example-block not executed in draft mode >>
```


## One harmonic {#One-harmonic}

The harmonic ansatz needs to be specified now – we expand `x` in a single frequency $\omega$.

```julia
add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)
```


```
<< @example-block not executed in draft mode >>
```


The object `diff_eq` now contains all the necessary information to convert the differential equation to the algebraic _harmonic equations_ (coupled polynomials in $U$ and $V$). 

```julia
harmonic_eq = get_harmonic_equations(diff_eq)
```


```
<< @example-block not executed in draft mode >>
```


The variables `u1` and `v1` were declared automatically to construct the harmonic ansatz. The _slow time_ variable `T` describes variation of the quadratures on timescales much slower than `ω`. For a steady state, all derivatives w.r.t `T` vanish, leaving only algebraic equations to be solved. 

We are ready to start plugging in numbers! Let us find steady states by solving `harmonic_eq` for numerical parameters. Homotopy continuation is especially suited to [solving over a range of parameter values](https://www.juliahomotopycontinuation.org/guides/parameter-homotopies/). Here we will solve over a range of driving frequencies `ω` – these are stored as `Pairs{Sym, Vector{Float}}`:

```julia
varied = ω => range(0.9, 1.2, 100); # range of parameter values
```


```
<< @example-block not executed in draft mode >>
```


The other parameters we be fixed – these are declared as `Pairs{Sym, Float}` pairs:

```julia
fixed = (α => 1., ω0 => 1.0, F => 0.01, γ => 0.01); # fixed parameters
```


```
<< @example-block not executed in draft mode >>
```


Now everything is ready to crank the handle. `get_steady_states` solves our `harmonic_eq` using the varied and fixed parameters:

```julia
result = get_steady_states(harmonic_eq, varied, fixed)
```


```
<< @example-block not executed in draft mode >>
```


The algorithm has found 3 solution branches in total (out of the [hypothetically admissible](https://en.wikipedia.org/wiki/B%C3%A9zout%27s_theorem) $3^{2} = 9$). All of these are real – and therefore physically observable – for at least some values of $\omega$. Only 2 branches are stable under infinitesimal perturbations. The &quot;Classes&quot; are boolean labels classifying each solution point, which may be used to select results for plotting.

We now want to visualize the results. Here we plot the solution amplitude, $\sqrt{U^2 + V^2}$ against the drive frequency $\omega$: 

```julia
plot(result, "sqrt(u1^2 + v1^2)")
```


```
<< @example-block not executed in draft mode >>
```


This is the expected [response curve](https://en.wikipedia.org/wiki/Duffing_equation#Frequency_response) for the Duffing equation.

## Using multiple harmonics {#Using-multiple-harmonics}

In the above section, we truncated the Fourier space to a single harmonic $\omega$ – the oscillator was assumed to only oscillate at the drive frequency. However, the Duffing oscillator can exhibit a rich spectrum of harmonics. We can obtain some intuition by treating $\alpha$ perturbatively in the equation of motion, i.e., by solving

$$\ddot{x}(t) + \gamma \dot{x}(t) + \omega_0^2 x(t) + \epsilon \alpha x(t)^3 = F \cos(\omega t)$$

for small $\epsilon$. To zeroth order, the response of the system is $x_0(t) = X_0 \cos(\omega t + \phi_0)$. Expanding $x(t) = x_0(t) + \epsilon x_1(t)$, we find that the perturbation $x_1(t)$ satisfies to first order

$$\ddot{x}_1(t) + \gamma \dot{x}_1(t) \left[ \omega_0^2 + \frac{3 \alpha X_0^2}{4} \right] x_1(t) = - \frac{\alpha X_0^3}{4} \cos(3 \omega t + 3 \phi_0) \,,$$

which gives a response of the form $x_1(t) = X_1 \cos(3 \omega t + \phi_1)$. Clearly, the oscillator now responds not only at frequency $\omega$, but also at $3 \omega$! This effect is known as [_high harmonic generation_](https://en.wikipedia.org/wiki/High_harmonic_generation) or more generally _frequency conversion_. By continuing the procedure to higher orders, we eventually obtain an infinity of harmonics present in the response. In general, there is no analytical solution to such problems.

We argued that frequency conversion takes place, to first order from $\omega$ to $3 \omega$. We can reflect this process by using a extended harmonic ansatz:

$$\begin{equation}
x(t) = U_1 \cos(\omega t) + V_1 \sin(\omega t) + U_2 \cos(3\omega t) + V_2 \sin(3\omega t) \,.
\end{equation}$$

Note that this is not a perturbative treatment! The harmonics $\omega$ and $3 \omega$ are on the same footing here. This is implemented as

```julia
add_harmonic!(diff_eq, x, [ω, 3ω]) # specify the two-harmonics ansatz
harmonic_eq = get_harmonic_equations(diff_eq)
```


```
<< @example-block not executed in draft mode >>
```


The variables `u1`, `v1` now encode `ω` and `u2`, `v2` encode `3ω`. We see this system is much harder to solve as we now have 4 harmonic variables, resulting in 4 coupled cubic equations. A maximum of $3^4 = 81$ solutions [may appear](https://en.wikipedia.org/wiki/B%C3%A9zout%27s_theorem)!

```julia
result = get_steady_states(harmonic_eq, varied, fixed)
plot(result, "sqrt(u1^2 + v1^2)")
```


```
<< @example-block not executed in draft mode >>
```


For the above parameters (where a perturbative treatment would have been reasonable), the principal response at $\omega$ looks rather similar, with a much smaller upconverted component appearing at $3 \omega$:

```julia
p1=plot(result, "sqrt(u1^2 + v1^2)", legend=false)
p2=plot(result, "sqrt(u2^2 + v2^2)")
plot(p1, p2)
```


```
<< @example-block not executed in draft mode >>
```


The non-perturbative nature of the ansatz allows us to capture some behaviour which is _not_ a mere extension of the usual single-harmonic Duffing response. Suppose we drive a strongly nonlinear resonator at frequency $\omega \cong \omega_0 / 3$. Such a drive is far out of resonance, however, the upconverted harmonic $3 \omega = \omega_0$ is not and may play an important role! Let us try this out:

```julia
fixed = (α => 10., ω0 => 3, F => 5, γ=>0.01)   # fixed parameters
varied = ω => range(0.9, 1.4, 100)           # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed)
```


```
<< @example-block not executed in draft mode >>
```


Although 9 branches were found in total, only 3 remain physical (real-valued). Let us visualise the amplitudes corresponding to the two harmonics, $\sqrt{U_1^2 + V_1^2}$ and $\sqrt{U_2^2 + V_2^2}$ :

```julia
p1 = plot(result, "sqrt(u1^2 + v1^2)", legend=false)
p2 = plot(result, "sqrt(u2^2 + v2^2)")
plot(p1, p2)
```


```
<< @example-block not executed in draft mode >>
```


The contributions of $\omega$ and $3\omega$ are now comparable and the system shows some fairly complex behaviour! This demonstrates how an exact solution within an extended Fourier subspace goes beyond a perturbative treatment.
