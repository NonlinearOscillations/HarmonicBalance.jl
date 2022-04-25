# [Introduction: linear response](@id linresp_ex)

In HarmonicBalance.jl, the [stability and linear response](@ref linresp_background) are treated using the [`LinearResponse`](@ref linresp_man) module. 

## Example: driven Duffing resonator

Here we calculate the white noise response of a simple nonlinear system. A set of reference results may be found in Huber et al. in [Phys. Rev. X 10, 021066 (2020)](https://doi.org/10.1103/PhysRevX.10.021066). The code is also available as a [Jupyter notebook](https://github.com/NonlinearOscillations/HarmonicBalance-notebooks). 

We start by [defining the Duffing oscillator](https://nonlinearoscillations.github.io/HarmonicBalance.jl/stable/examples/single_Duffing/#The-code)
```julia
using HarmonicBalance
@variables α, ωF, ω0, F, t, q(t), γ, Γ; # declare constant variables and a function q(t)

# define ODE
diff_eq = DifferentialEquation(d(q,t,2) + 2*Γ*d(q,t) + ω0^2*q + γ*q^3 ~ F*cos(ωF*t), q)

# specify the ansatz q = u(T) cos(ωF*t) + v(T) sin(ωF*t)
add_harmonic!(diff_eq, q, ωF) 

# implement ansatz to get harmonic equations
harmonic_eq = get_harmonic_equations(diff_eq)
```

### Linear regime

When driven weakly, the Duffing resonator behaves quasi-linearly, i.e, its response to noise is independent of the applied drive. We see that for weak driving, $F = 10^{-6}$, the amplitude is a Lorentzian. 
```julia
fixed = (α => 1, ω0 => 1.0, γ => 1E-2, F => 1E-6)   # fixed parameters
swept = ω => LinRange(0.9, 1.1, 100)           # range of parameter values
solutions = get_steady_states(harmonic_eq, swept, fixed)

plot_1D_solutions(solutions, x="ω", y="sqrt(u1^2 + v1^2)")
```

```@raw html
<img style="display: block; margin: 0 auto; padding-bottom: 20px;" src="../../assets/linear_response/Duffing_quasilin_amp.png" width="450" alignment="left" \>
``` ⠀
The linear response is obtained with [`plot_jacobian_spectrum`](@ref HarmonicBalance.LinearResponse.plot_jacobian_spectrum), note that the branch number and the variable (here `x`) must be specified. The response has a peak at $\omega_0$, irrespective of $\omega$ :
```julia
HarmonicBalance.LinearResponse.plot_jacobian_spectrum(solutions, x, 
    Ω_range=LinRange(0.9, 1.1, 300), branch=1, logscale=true)
```
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/linear_response/Duffing_quasilin_noise.png" width="450" alignment="left" \>
``` ⠀
Note the slight "bending" of the noise peak with $\omega$ - this is given by the failure of the Jacobian to capture response far-detuned from the drive frequency. More on this point will follow in a future example.

### Nonlinear regime

For strong driving, matters get more complicated. Let us now use a drive $F = 10^{-2}$ :

```julia
fixed = (α => 1, ω0 => 1.0, γ => 1E-2, F => 1E-2)   # fixed parameters
swept = ω => LinRange(0.9, 1.1, 100)           # range of parameter values
solutions = get_steady_states(harmonic_eq, swept, fixed)

plot_1D_solutions(solutions, x="ω", y="sqrt(u1^2 + v1^2)");
```
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/linear_response/Duffing_nonlin_amp.png" width="450" alignment="left" \>
``` ⠀

The amplitude is the well-known Duffing curve. Let's see the linear response of the two stable branches, 1 and 2.
```julia
LinearResponse.plot_jacobian_spectrum(solutions, x, 
    Ω_range=LinRange(0.9,1.1,300), branch=1, logscale=true);

LinearResponse.plot_jacobian_spectrum(solutions, x, 
    Ω_range=LinRange(0.9,1.1,300), branch=2, logscale=true);
```

```@raw html
<img style="display: block; padding-bottom: 20px;" src="../../assets/linear_response/Duffing_nonlin_noise12.png" width="900" alignment="center" \>
``` ⠀
In branch 1 the linear response to white noise shows _more than one peak_. This is a distinctly nonlinear phenomenon, manifesting primarily at large amplitudes. Branch 2 is again quasi-linear, which stems from its low amplitude.

Following [Huber et al.](https://doi.org/10.1103/PhysRevX.10.021066), we may also fix $\omega = \omega_0$ and plot the linear response as a function of $F$. The response turns out to be single-valued over a large range of driving strengths. Using a log scale for the x-axis:

```julia
fixed = (α => 1., ω0 => 1.0, γ => 1E-2, ω => 1)   # fixed parameters
swept = F => 10 .^ LinRange(-6, -1, 200)           # range of parameter values
solutions = get_steady_states(harmonic_eq, swept, fixed)

plot_1D_solutions(solutions, x="F", y="sqrt(u1^2 + v1^2)");
HarmonicBalance.xscale("log") # use log scale on x

LinearResponse.plot_jacobian_spectrum(solutions, x, 
    Ω_range=LinRange(0.9,1.1,300), branch=1, logscale=true);
HarmonicBalance.xscale("log")
```
```@raw html
<img style="display: block; padding-bottom: 20px;" src="../../assets/linear_response/Duffing_nonlin_forcex.png" width="800" alignment="center" \>
``` ⠀
We see that for low $F$, quasi-linear behaviour with a single Lorentzian response occurs, while for larger $F$, two peaks form in the noise response. The two peaks are strongly unequal in magnitude, which is an example of internal squeezing.



