


# Parametric Pumping via Three-Wave Mixing {#Parametric-Pumping-via-Three-Wave-Mixing}

```julia
using HarmonicBalance, Plots
using Plots.Measures
using Random
```


```
<< @example-block not executed in draft mode >>
```


## System

```julia
@variables β α ω ω0 F γ t x(t) # declare constant variables and a function x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + β * x^2 + α * x^3 + γ * d(x, t) ~ F * cos(ω * t), x
)
add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)
```


```
<< @example-block not executed in draft mode >>
```


## 1st order Krylov expansion {#1st-order-Krylov-expansion}

```julia
harmonic_eq = get_krylov_equations(diff_eq; order=1)
harmonic_eq.equations
```


```
<< @example-block not executed in draft mode >>
```


If we both have quadratic and cubic nonlineariy, we observe the normal duffing oscillator response.

```julia
varied = (ω => range(0.99, 1.1, 200)) # range of parameter values
fixed = (α => 1.0, β => 1.0, ω0 => 1.0, γ => 0.005, F => 0.0025) # fixed parameters

result = get_steady_states(harmonic_eq, varied, fixed)
plot(result; y="u1^2+v1^2")
```


```
<< @example-block not executed in draft mode >>
```


If we set the cubic nonlinearity to zero, we recover the driven damped harmonic oscillator. Indeed, thefirst order the quadratic nonlinearity has no affect on the system.

```julia
varied = (ω => range(0.99, 1.1, 100))
fixed = (α => 0.0, β => 1.0, ω0 => 1.0, γ => 0.005, F => 0.0025)

result = get_steady_states(harmonic_eq, varied, fixed)
plot(result; y="u1^2+v1^2")
```


```
<< @example-block not executed in draft mode >>
```


## 2nd order Krylov expansion {#2nd-order-Krylov-expansion}

The quadratic nonlinearity $\beta$ together with the drive at 2ω gives the effective parametric drive $\lambda_\mathrm{eff}=\frac{2F_1\beta}{3m\omega^2}$. But the cubic nonlinearity $\alpha$ is still needed to get the period doubling bifurcation through $\lambda_\mathrm{eff}$.

```julia
@variables β α ω ω0 F γ t x(t)
diff_eq = DifferentialEquation(
    d(x, t, 2) + ω0^2 * x + β * x^2 + α * x^3 + γ * d(x, t) ~ F * cos(2ω * t), x
)

add_harmonic!(diff_eq, x, ω)
harmonic_eq2 = get_krylov_equations(diff_eq; order=2)
```


```
<< @example-block not executed in draft mode >>
```


```julia
varied = (ω => range(0.4, 1.1, 500))
fixed = (α => 1.0, β => 2.0, ω0 => 1.0, γ => 0.001, F => 0.005)

result = get_steady_states(harmonic_eq2, varied, fixed)
plot(result; y="v1")
```


```
<< @example-block not executed in draft mode >>
```


```julia
varied = (ω => range(0.4, 0.6, 100), F => range(1e-6, 0.01, 50))
fixed = (α => 1.0, β => 2.0, ω0 => 1.0, γ => 0.01)

method = TotalDegree()
result = get_steady_states(harmonic_eq2, method, varied, fixed)
plot_phase_diagram(result; class="stable")
```


```
<< @example-block not executed in draft mode >>
```



---


_This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl)._
