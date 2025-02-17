# [Linear response](@id linresp_ex)

In HarmonicBalance.jl, the [stability and linear response](@ref linresp_background) are treated using the [`LinearResponse`](@ref linresp_man) module.

Here we calculate the white noise response of a simple nonlinear system. A set of reference results may be found in Huber et al. in [Phys. Rev. X 10, 021066 (2020)](https://doi.org/10.1103/PhysRevX.10.021066).
We start by defining the [Duffing oscillator](@ref Duffing)

```@example linresp
using HarmonicBalance
using Plots.Measures: mm
@variables α, ω, ω0, F, γ, t, x(t); # declare constant variables and a function x(t)

# define ODE
diff_eq = DifferentialEquation(d(x,t,2) + ω0*x + α*x^3 + γ*d(x,t) ~ F*cos(ω*t), x)

# specify the ansatz x = u(T) cos(ω*t) + v(T) sin(ω*t)
add_harmonic!(diff_eq, x, ω) 

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq)
```

## Linear regime

When driven weakly, the Duffing resonator behaves quasi-linearly, i.e, its response to noise is independent of the applied drive. We see that for weak driving, $F = 10^{-4}$, the amplitude is a Lorentzian.

```@example linresp
fixed = (α => 1, ω0 => 1.0, γ => 0.005, F => 0.0001)   # fixed parameters
varied = ω => range(0.95, 1.05, 100)           # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed)

using Plots
plot(result, "sqrt(u1^2 + v1^2)")
```

To find the fluctuation on the top of the steady state one often employs a [Bogoliubov-de Gennes analyses](https://en.wikipedia.org/wiki/Linear_dynamical_system). Here, we compute the eigenvalues $\lambda_k$ of the Jacobian matrix at the steady state. The imaginary part of the eigenvalues gives characteristic frequencies of the "quasi-particle excitations". The real part gives the lifetime of these excitations.

The compute the eigenvalues of a specific branch, we can use the corresponding function:

```@example linresp
eigvalues= eigenvalues(result, 1)
```

Using the [PlotsExt.jl](@ref plotting) extension, one can quickly compute and plot the eigenvalues as follows

```@example linresp
plot(
    plot_eigenvalues(result, 1),
    plot_eigenvalues(result, 1, type=:real, ylims=(-0.003, 0)),
)
```

We find a single pair of complex conjugate eigenvalues linearly changing with the driving frequency. Both real parts are negative, indicating stability.

As discussed in [background section on linear response](@ref linresp_background), the excitation manifest itself as a lorentenzian peak in a Power Spectral Density (PSD) measurement. The PSD can be plotted using [`plot_linear_response`](@ref linresp_man):

```@example linresp
plot_linear_response(result, x, 1, Ω_range=range(0.95, 1.05, 300), logscale=true)
```

The response has a peak at $\omega_0$, irrespective of the driving frequency $\omega$. Indeed, the eigenvalues shown before where plotted in the rotating frame at the frequency of the drive $\omega$. Hence, the imaginary part of eigenvalues shows the frequency (energy) needed to excite the system at it natural frequency (The frequency its want to be excited at.)

Note the slight "bending" of the noise peak with $\omega$ - this is given by the failure of the first-order calculation of the jacobian to capture response far-detuned from the drive frequency. One can correct this by using higher-order derivatives of the `Differentialequation` object in the jacobian calculation. For more details on this see the [thesis](https://www.doi.org/10.3929/ethz-b-000589190). We can use this corrections by setting the `order` argument in the `plot_linear_response` function:

```@example linresp
plot_linear_response(result, x, 1, Ω_range=range(0.95, 1.05, 300), logscale=true, order=2)
```

To compute the matrix without plotting you can use the functions specified at the [linear respinse manual](@ref linresp_man).

## Nonlinear regime

For strong driving, matters get more complicated. Let us now use a drive $F = 2*10^{-3}$ :

```@example linresp
fixed = (α => 1, ω0 => 1.0, γ => 0.005, F => 0.002)   # fixed parameters
varied = ω => range(0.95, 1.05, 100)           # range of parameter values
result = get_steady_states(harmonic_eq, varied, fixed)

plot(result, x="ω", y="sqrt(u1^2 + v1^2)");
```

The amplitude is the well-known Duffing curve. Let's look at the eigenvalues of the two stable branches, 1 and 2.

```@example linresp
plot(
    plot_eigenvalues(result, 1),
    plot_eigenvalues(result, 1, type=:real, ylims=(-0.003, 0)),
    plot_eigenvalues(result, 2),
    plot_eigenvalues(result, 2, type=:real, ylims=(-0.003, 0)),
)
```

Again every branch gives a single pair of complex conjugate eigenvalues. However, for branch 1, the characteristic frequencies due not change linearly with the driving frequency around $\omega=\omega_0$. This is a sign of steady state becoming nonlinear at large amplitudes.

The same can be seen in the PSD:

```@example linresp
plot(
  plot_linear_response(result, x, 1, Ω_range=range(0.95,1.1,300), logscale=true),
  plot_linear_response(result, x, 2, Ω_range=range(0.9,1.1,300), logscale=true),
    size=(600, 250), margin=3mm
)
```

In branch 1 the linear response to white noise shows _more than one peak_. This is a distinctly nonlinear phenomenon, indicative of the squeezing of the steady state. Branch 2 is again quasi-linear, which stems from its low amplitude.

We can compute the squeezing of the steady states by using the corresponding eigenvectors of the eigenvalus. Indeed, defining (TODO add reference)

```@example linresp
function symplectic(v)
    2 * (real(v[1]) * imag(v[2]) - imag(v[1]) * real(v[2]))
end
function squeeze(v)
    symp = symplectic(v)
    ((1 - symp) / (1 + symp))^sign(symp)
end
```

We can compute the squeezing of the steady states as follows:

```@example linresp
eigvecs = eigenvectors(result, 1)
squeezed = [squeeze.(eachcol(mat))[1] for mat in eigvecs]
plot(range(0.95, 1.05, 100), squeezed, label="Squeezing of branch 1")
```

Following [Huber et al.](https://doi.org/10.1103/PhysRevX.10.021066), we may also fix $\omega = \omega_0$ and plot the linear response as a function of $F$. The response turns out to be single-valued over a large range of driving strengths. Using a log scale for the x-axis:

```@example linresp
fixed = (α => 1., ω0 => 1.0, γ => 1e-2, ω => 1)   # fixed parameters
swept = F => 10 .^ range(-6, -1, 200)           # range of parameter values
result = get_steady_states(harmonic_eq, swept, fixed)

plot(
  plot(result, "sqrt(u1^2 + v1^2)", xscale=:log),
  plot_linear_response(result, x, 1, Ω_range=range(0.9,1.1,300), logscale=true, xscale=:log),
  size=(600, 250), margin=3mm
)
```

We see that for low $F$, quasi-linear behaviour with a single Lorentzian response occurs, while for larger $F$, two peaks form in the noise response. The two peaks are strongly unequal in magnitude, which is an example of internal squeezing (See supplemental material of [Huber et al.](https://doi.org/10.1103/PhysRevX.10.021066)).
