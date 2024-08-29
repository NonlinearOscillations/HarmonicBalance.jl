# [Introduction: linear response](@id linresp_ex)

In HarmonicBalance.jl, the [stability and linear response](@ref linresp_background) are treated using the [`LinearResponse`](@ref linresp_man) module. 

Here we calculate the white noise response of a simple nonlinear system. A set of reference results may be found in Huber et al. in [Phys. Rev. X 10, 021066 (2020)](https://doi.org/10.1103/PhysRevX.10.021066).
We start by [defining the Duffing oscillator](https://nonlinearoscillations.github.io/HarmonicBalance.jl/stable/examples/simple_Duffing/#The-code)
```@example linresp
using HarmonicBalance, Plots
using Plots.Measures: mm
@variables α, ω, ω0, F, γ, t, x(t); # declare constant variables and a function x(t)

# define ODE
diff_eq = DifferentialEquation(d(x,t,2) + ω0*x + α*x^3 + γ*d(x,t) ~ F*cos(ω*t), x)

# specify the ansatz x = u(T) cos(ω*t) + v(T) sin(ω*t)
add_harmonic!(diff_eq, x, ω) 

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq)
```

### Linear regime

When driven weakly, the Duffing resonator behaves quasi-linearly, i.e, its response to noise is independent of the applied drive. We see that for weak driving, $F = 10^{-6}$, the amplitude is a Lorentzian. 
```@example linresp
fixed = (α => 1, ω0 => 1.0, γ => 1e-2, F => 1e-6)   # fixed parameters
varied = ω => range(0.9, 1.1, 100)           # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed)

plot(result, "sqrt(u1^2 + v1^2)")
```
The linear response is obtained with [`plot_linear_response`](@ref linresp_man), note that the branch number and the variable (here `x`) must be specified. The response has a peak at $\omega_0$, irrespective of $\omega$ :
```@example linresp
plot_linear_response(result, x, Ω_range=range(0.9,1.1,300), branch=1, logscale=true)
```
Note the slight "bending" of the noise peak with $\omega$ - this is given by the failure of the first-order calculation to capture response far-detuned from the drive frequency. More on this point will follow in a future example.

### Nonlinear regime

For strong driving, matters get more complicated. Let us now use a drive $F = 10^{-2}$ :

```@example linresp
fixed = (α => 1, ω0 => 1.0, γ => 1e-2, F => 1e-2)   # fixed parameters
swept = ω => range(0.9, 1.1, 100)           # range of parameter values
result = get_steady_states(harmonic_eq, swept, fixed)

plot(result, x="ω", y="sqrt(u1^2 + v1^2)");
```

The amplitude is the well-known Duffing curve. Let's see the linear response of the two stable branches, 1 and 2.
```@example linresp
plot(
  plot_linear_response(result, x, branch=1, Ω_range=range(0.9,1.1,300), logscale=true),
  plot_linear_response(result, x, branch=2, Ω_range=range(0.9,1.1,300), logscale=true),
    size=(600, 250), margin=3mm
)
```
In branch 1 the linear response to white noise shows _more than one peak_. This is a distinctly nonlinear phenomenon, manifesting primarily at large amplitudes. Branch 2 is again quasi-linear, which stems from its low amplitude.

Following [Huber et al.](https://doi.org/10.1103/PhysRevX.10.021066), we may also fix $\omega = \omega_0$ and plot the linear response as a function of $F$. The response turns out to be single-valued over a large range of driving strengths. Using a log scale for the x-axis:

```@example linresp
fixed = (α => 1., ω0 => 1.0, γ => 1e-2, ω => 1)   # fixed parameters
swept = F => 10 .^ range(-6, -1, 200)           # range of parameter values
result = get_steady_states(harmonic_eq, swept, fixed)

plot(
  plot(result, "sqrt(u1^2 + v1^2)", xscale=:log),
  plot_linear_response(result, x, branch=1, Ω_range=range(0.9,1.1,300), logscale=true, xscale=:log),
  size=(600, 250), margin=3mm
)
```
We see that for low $F$, quasi-linear behaviour with a single Lorentzian response occurs, while for larger $F$, two peaks form in the noise response. The two peaks are strongly unequal in magnitude, which is an example of internal squeezing.



