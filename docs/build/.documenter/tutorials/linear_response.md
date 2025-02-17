
# Linear response {#linresp_ex}

In HarmonicBalance.jl, the [stability and linear response](/background/stability_response#linresp_background) are treated using the [`LinearResponse`](/manual/linear_response#linresp_man) module. 

Here we calculate the white noise response of a simple nonlinear system. A set of reference results may be found in Huber et al. in [Phys. Rev. X 10, 021066 (2020)](https://doi.org/10.1103/PhysRevX.10.021066). We start by defining the [Duffing oscillator](/tutorials/steady_states#Duffing)

```julia
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


```
<< @example-block not executed in draft mode >>
```


### Linear regime {#Linear-regime}

When driven weakly, the Duffing resonator behaves quasi-linearly, i.e, its response to noise is independent of the applied drive. We see that for weak driving, $F = 10^{-4}$, the amplitude is a Lorentzian. 

```julia
fixed = (α => 1, ω0 => 1.0, γ => 0.005, F => 0.0001)   # fixed parameters
varied = ω => range(0.95, 1.05, 100)           # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed)

plot(result, "sqrt(u1^2 + v1^2)")
```


```
<< @example-block not executed in draft mode >>
```


To find the fluctuation on the top of the steady state one often employs a [Bogoliubov-de Gennes analyses](https://en.wikipedia.org/wiki/Linear_dynamical_system). Here, we compute the eigenvalues $\lambda_k$ of the Jacobian matrix at the steady state. The imaginary part of the eigenvalues gives characteristic frequencies of the &quot;quasi-particle excitations&quot;. The real part gives the lifetime of these excitations.

One can plot the eigenvalues as follows

```julia
plot(
    plot_eigenvalues(result, branch=1),
    plot_eigenvalues(result, branch=1, type=:real, ylims=(-0.003, 0)),
)
```


```
<< @example-block not executed in draft mode >>
```


We find a single pair of complex conjugate eigenvalues linearly changing with the driving frequency. Both real parts are negative, indicating stability.

As discussed in [background section on linear response](/background/stability_response#linresp_background), the excitation manifest itself as a lorentenzian peak in a power spectral density (PSD) measurement. The PSD can be plotted using [`plot_linear_response`](/manual/linear_response#linresp_man):

```julia
plot_linear_response(result, x, Ω_range=range(0.95, 1.05, 300), branch=1, logscale=true)
```


```
<< @example-block not executed in draft mode >>
```


The response has a peak at $\omega_0$, irrespective of the driving frequency $\omega$. Indeed, the eigenvalues shown before where plotted in the rotating frame at the frequency of the drive $\omega$. Hence, the imaginary part of eigenvalues shows the frequency (energy) needed to excite the system at it natural frequency (The frequency its want to be excited at.)

Note the slight &quot;bending&quot; of the noise peak with $\omega$ - this is given by the failure of the first-order calculation to capture response far-detuned from the drive frequency.

### Nonlinear regime {#Nonlinear-regime}

For strong driving, matters get more complicated. Let us now use a drive $F = 2*10^{-3}$ :

```julia
fixed = (α => 1, ω0 => 1.0, γ => 0.005, F => 0.002)   # fixed parameters
varied = ω => range(0.95, 1.05, 100)           # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed)

plot(result, x="ω", y="sqrt(u1^2 + v1^2)");
```


```
<< @example-block not executed in draft mode >>
```


The amplitude is the well-known Duffing curve. Let&#39;s look at the eigenvalues of the two stable branches, 1 and 2.

```julia
plot(
    plot_eigenvalues(result, branch=1),
    plot_eigenvalues(result, branch=1, type=:real, ylims=(-0.003, 0)),
    plot_eigenvalues(result, branch=2),
    plot_eigenvalues(result, branch=2, type=:real, ylims=(-0.003, 0)),
)
```


```
<< @example-block not executed in draft mode >>
```


Again every branch gives a single pair of complex conjugate eigenvalues. However, for branch 1, the characteristic frequencies due not change linearly with the driving frequency around $\omega=\omega_0$. This is a sign of steady state becoming nonlinear at large amplitudes.

The same can be seen in the PSD:

```julia
plot(
  plot_linear_response(result, x, branch=1, Ω_range=range(0.95,1.1,300), logscale=true),
  plot_linear_response(result, x, branch=2, Ω_range=range(0.9,1.1,300), logscale=true),
    size=(600, 250), margin=3mm
)
```


```
<< @example-block not executed in draft mode >>
```


In branch 1 the linear response to white noise shows _more than one peak_. This is a distinctly nonlinear phenomenon, indicative if the squeezing of the steady state. Branch 2 is again quasi-linear, which stems from its low amplitude.

Following [Huber et al.](https://doi.org/10.1103/PhysRevX.10.021066), we may also fix $\omega = \omega_0$ and plot the linear response as a function of $F$. The response turns out to be single-valued over a large range of driving strengths. Using a log scale for the x-axis:

```julia
fixed = (α => 1., ω0 => 1.0, γ => 1e-2, ω => 1)   # fixed parameters
swept = F => 10 .^ range(-6, -1, 200)           # range of parameter values
result = get_steady_states(harmonic_eq, swept, fixed)

plot(
  plot(result, "sqrt(u1^2 + v1^2)", xscale=:log),
  plot_linear_response(result, x, branch=1, Ω_range=range(0.9,1.1,300), logscale=true, xscale=:log),
  size=(600, 250), margin=3mm
)
```


```
<< @example-block not executed in draft mode >>
```


We see that for low $F$, quasi-linear behaviour with a single Lorentzian response occurs, while for larger $F$, two peaks form in the noise response. The two peaks are strongly unequal in magnitude, which is an example of internal squeezing (See supplemental material of [Huber et al.](https://doi.org/10.1103/PhysRevX.10.021066)).
