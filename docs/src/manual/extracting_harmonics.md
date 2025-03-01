# Extracting harmonic equations

## [Harmonic Balance method](@id Harmonic_Balance)

Once a `DifferentialEquation` is defined and its harmonics specified, one can extract the harmonic equations using `get_harmonic_equations`, which itself is composed of the subroutines `harmonic_ansatz`, `slow_flow`, `fourier_transform!` and `drop_powers`.

The harmonic equations use an additional time variable specified as `slow_time` in `get_harmonic_equations`. This is essentially a label distinguishing the time dependence of the harmonic variables (expected to be slow)
from that of the oscillating terms (expected to be fast). When the equations are Fourier-transformed to remove oscillating terms, `slow_time` is treated as a constant. Such an approach is exact when looking for steady states.

```@docs; canonical=false
get_harmonic_equations
```

## HarmonicVariable and HarmonicEquation types

The equations governing the harmonics are stored using the two following structs. When going from the original to the harmonic equations, the harmonic ansatz $x_i(t) = \sum_{j=1}^M u_{i,j}  (T)  \cos(\omega_{i,j} t)+ v_{i,j}(T) \sin(\omega_{i,j} t)$ is used. Internally, each pair $(u_{i,j}, v_{i,j})$ is stored as a `HarmonicVariable`. This includes the identification of $\omega_{i,j}$ and $x_i(t)$, which is needed to later reconstruct $x_i(t)$.

```@docs; canonical=false
HarmonicVariable
```

When the full set of equations of motion is expanded using the harmonic ansatz, the result is stored as a `HarmonicEquation`. For an initial equation of motion consisting of $M$ variables, each expanded in $N$ harmonics, the resulting `HarmonicEquation` holds $2NM$ equations of $2NM$ variables. Each symbol not corresponding to a variable is identified as a parameter.

A `HarmonicEquation` can be either parsed into a steady-state [`HarmonicBalance.Problem`](@ref) or solved using a dynamical ODE solver.

```@docs; canonical=false
HarmonicEquation
```

## [Krylov-Bogoliubov Averaging Method](@id Krylov-Bogoliubov)

The Krylov-Bogoliubov averaging method is an alternative high-frequency expansion technique used to analyze dynamical systems. Unlike the [Harmonic Balance method](https://en.wikipedia.org/wiki/Harmonic_balance), which is detailed in the [background section](@ref intro_hb), the Krylov-Bogoliubov method excels in computing higher orders in $1/\omega$, enabling the capture of faster dynamics within a system.

```@docs; canonical=false
get_krylov_equations
```

### Purpose and Advantages

The primary advantage of the Krylov-Bogoliubov method lies in its ability to delve deeper into high-frequency components, allowing a more comprehensive understanding of fast dynamical behaviors. By leveraging this technique, one can obtain higher-order approximations that shed light on intricate system dynamics.

However, it's essential to note a limitation: this method cannot handle multiple harmonics within a single variable, unlike some other high-frequency expansion methods.

For further information and a detailed understanding of this method, refer to [Krylov-Bogoliubov averaging method on Wikipedia](https://en.wikipedia.org/wiki/Krylov%E2%80%93Bogoliubov_averaging_method).
