
# Classifying solutions {#classifying}

Given that you obtained some steady states for a parameter sweep of a specific model it can be useful to classify these solution. Let us consider a simple pametric oscillator 

```julia
using HarmonicBalance

@variables ω₀ γ λ α ω t x(t)

natural_equation = d(d(x, t), t) + γ * d(x, t) + (ω₀^2 - λ * cos(2 * ω * t)) * x + α * x^3
diff_eq = DifferentialEquation(natural_equation, x)

add_harmonic!(diff_eq, x, ω);

harmonic_eq = get_harmonic_equations(diff_eq)
```


```
<< @example-block not executed in draft mode >>
```


We perform a 2d sweep in the driving frequency $\omega$ and driving strength $\lambda$:

```julia
fixed = (ω₀ => 1.0, γ => 0.002, α => 1.0)
varied = (ω => range(0.99, 1.01, 100), λ => range(1e-6, 0.03, 100))

result_2D = get_steady_states(harmonic_eq, varied, fixed)
```


```
<< @example-block not executed in draft mode >>
```


By default the steady states of the system are classified by four different catogaries:
- `physical`: Solutions that are physical, i.e., all variables are purely real.
  
- `stable`: Solutions that are stable, i.e., all eigenvalues of the Jacobian have negative real parts.
  
- `Hopf`: Solutions that are physical and have exactly two Jacobian eigenvalues with positive real parts, which are complex conjugates of each other. The class can help to identify regions where a limit cycle is present due to a [Hopf bifurcation](https://en.wikipedia.org/wiki/Hopf_bifurcation). See also the tutorial on [limit cycles](/tutorials/limit_cycles#limit_cycles).
  
- `binary_labels`: each region in the parameter sweep receives an identifier based on its permutation of stable branches. This allows to distinguish between different phases, which may have the same number of stable solutions.
  

We can plot the number of stable solutions, giving the phase diagram

```julia
plot_phase_diagram(result_2D, class="stable")
```


```
<< @example-block not executed in draft mode >>
```


If we plot the a cut at $\lambda=0.01$, we see that in the blue region only one stable solution exists with zero amplitude:

```julia
plot(result_2D, y="√(u1^2+v1^2)", cut=λ => 0.01, class="stable") |> display
```


```
<< @example-block not executed in draft mode >>
```


```julia
get_single_solution(result_2D; branch=1, index=(1, 1))
```


```
<< @example-block not executed in draft mode >>
```


This solution becomes stable again outside the green lobe. Also called Mathieu lobe. Indeed, we can classify the zero amplitude solution by adding an extra category as a class:

```julia
classify_solutions!(result_2D, "sqrt(u1^2 + v1^2) < 0.001", "zero")
result_2D
```


```
<< @example-block not executed in draft mode >>
```


We can visualize the zero amplitude solution:

```julia
plot_phase_diagram(result_2D, class=["zero", "stable"])
```


```
<< @example-block not executed in draft mode >>
```


This shows that inside the Mathieu lobe the zero amplitude solution becomes unstable due to the parametric drive being resonant with the oscillator.

We can also visualize the equi-amplitude curves of the solutions:

```julia
classify_solutions!(result_2D, "sqrt(u1^2 + v1^2) > 0.12", "large amplitude")
plot_phase_diagram(result_2D, class=["large amplitude", "stable"])
```


```
<< @example-block not executed in draft mode >>
```

