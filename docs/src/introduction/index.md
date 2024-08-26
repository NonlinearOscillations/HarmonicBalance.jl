Let's find the steady states of a driven Duffing oscillator with nonlinear damping, its equation of motion is:
```math
\begin{equation} \label{eq:duffing}
\underbrace{\ddot{x}(t) + \gamma \dot{x}(t) + \omega_0^2 x(t)}_{\text{damped harmonic oscillator}} + \underbrace{\alpha x(t)^3}_{\text{Duffing coefficient}} = \underbrace{F \cos(\omega t)}_{\text{periodic drive}}
\end{equation}
```

```julia
using HarmonicBalance
@variables α ω ω0 F t η x(t) # declare constant variables and a function x(t)
diff_eq = DifferentialEquation(d(x,t,2) + ω0^2*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t), x)
add_harmonic!(diff_eq, x, ω) # specify the ansatz x = u(T) cos(ωt) + v(T) sin(ωt)

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq)

fixed = (α => 1.0, ω0 => 1.0, F => 0.01, η => 0.1)   # fixed parameters
varied = ω => range(0.9, 1.2, 100)           # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed)
```
The results are `show`n:
```
A steady state result for 100 parameter points

Solution branches:   3
   of which real:    3
   of which stable:  2

Classes: stable, physical, Hopf, binary_labels
```
```julia
plot(result, "sqrt(u1^2 + v1^2)")
```
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/simple_Duffing/response_single.png" alignment="center" \>
``` ⠀
