# Extracting harmonic equations

## Harmonic Balance method
Once a `DifferentialEquation` is defined and its harmonics specified, one can extract the harmonic equations using `get_harmonic_equations`, which itself is composed of the subroutines `harmonic_ansatz`, `slow_flow`, `fourier_transform!` and `drop_powers`. 

The harmonic equations use an additional time variable specified as `slow_time` in `get_harmonic_equations`. This is essentially a label distinguishing the time dependence of the harmonic variables (expected to be slow)
from that of the oscillating terms (expeted to be fast). When the equations are Fourier-transformed to remove oscillating terms, `slow_time` is treated as a constant. Such an approach is exact when looking for steady states. 

```@docs
get_harmonic_equations
HarmonicBalance.harmonic_ansatz
HarmonicBalance.slow_flow
HarmonicBalance.fourier_transform
HarmonicBalance.drop_powers
```

## HarmonicVariable and HarmonicEquation types

The equations governing the harmonics are stored using the two following structs. When going from the original to the harmonic equations, the harmonic ansatz $x_i(t) = \sum_{j=1}^M u_{i,j}  (T)  \cos(\omega_{i,j} t)+ v_{i,j}(T) \sin(\omega_{i,j} t)$ is used. Internally, each pair $(u_{i,j}, v_{i,j})$ is stored as a `HarmonicVariable`. This includes the identification of $\omega_{i,j}$ and $x_i(t)$, which is needed to later reconstruct $x_i(t)$.

```@docs
HarmonicVariable
```

When the full set of equations of motion is expanded using the harmonic ansatz, the result is stored as a `HarmonicEquation`. For an initial equation of motion consisting of $M$ variables, each expanded in $N$ harmonics, the resulting `HarmonicEquation` holds $2NM$ equations of $2NM$ variables. Each symbol not corresponding to a variable is identified as a parameter. 

A `HarmonicEquation` can be either parsed into a steady-state `Problem` or solved using a dynamical ODE solver.

```@docs
HarmonicEquation
```