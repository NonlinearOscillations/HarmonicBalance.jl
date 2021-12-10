# Extracting harmonic equations

Once a `DifferentialEquation` is defined and its harmonics specified, one can extract the harmonic equations using `get_harmonic_equations`, which itself is composed of the subroutines `harmonic_ansatz`, `slow_flow`, `fourier_transform!` and `drop_powers`. 


```@docs
get_harmonic_equations
HarmonicBalance.harmonic_ansatz
HarmonicBalance.slow_flow
HarmonicBalance.fourier_transform
HarmonicBalance.drop_powers
```

## HarmonicVariable and HarmonicEquation types

The equations governing the harmonics are stored using the two following structs. When going from the original to the harmonic equations, the harmonic ansatz $x_i(t) = \sum_{j=1}^M u_i^{(j)}  (T)  \cos(\omega_j^{(i)} t)+ v_i^{(j)} (T) \sin(\omega_j^{(i)} t)$ is used. Internally, each pair $(u_i^{(j)}, v_i^{(j)})$ is stored as a `HarmonicVariable`. This includes the identification of $\omega_j^{(i)}$ and $x_i(t)$, which is needed to later reconstruct the solutions in the non-rotating frame.

```@docs
HarmonicVariable
```

When the full set of equations of motion is expanded using the harmonic ansatz, the result is stored as a `HarmonicEquation`. For an initial equation of motion consisting of $M$ variables, each expanded in $N$ harmonics, the resulting `HarmonicEquation` holds $2NM$ equations of $2NM$ variables. Each symbol not corresponding to a variable is identified as a parameter. 

A `HarmonicEquation` can be either parsed into a steady-state `Problem` or solved using a dynamical ODE solver.

```@docs
HarmonicEquation
```