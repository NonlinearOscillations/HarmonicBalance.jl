# [Krylov-Bogoliubov Averaging Method](@id Krylov-Bogoliubov)

The Krylov-Bogoliubov averaging method is an alternative high-frequency expansion technique used to analyze dynamical systems. Unlike the [Harmonic Balance method](https://en.wikipedia.org/wiki/Harmonic_balance), which is detailed in the [background section](@ref intro_hb), the Krylov-Bogoliubov method excels in computing higher orders in $1/\omega$, enabling the capture of faster dynamics within a system.

## Purpose and Advantages

The primary advantage of the Krylov-Bogoliubov method lies in its ability to delve deeper into high-frequency components, allowing a more comprehensive understanding of fast dynamical behaviors. By leveraging this technique, one can obtain higher-order approximations that shed light on intricate system dynamics.

However, it's essential to note a limitation: this method cannot handle multiple harmonics within a single variable, unlike some other high-frequency expansion methods.

## Usage

To compute the Krylov-Bogoliubov averaging method within your system, utilize the function `get_krylov_equations`. This function is designed specifically to implement the methodology and derive the equations necessary to analyze the system dynamics using this technique.

### Function Reference

```@docs; canonical=false
get_krylov_equations
```

For further information and a detailed understanding of this method, refer to [Krylov-Bogoliubov averaging method on Wikipedia](https://en.wikipedia.org/wiki/Krylov%E2%80%93Bogoliubov_averaging_method).
