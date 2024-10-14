# Entering equations of motion 

The struct `DifferentialEquation` is the primary input method; it holds an ODE or a coupled system of ODEs composed of terms with harmonic time-dependence
The dependent variables are specified during input, any other symbols
are identified as parameters. Information on which variable is to be expanded in which harmonic is specified using `add_harmonic!`. 

`DifferentialEquation.equations` stores a dictionary assigning variables to equations. This information is necessary because the harmonics belonging to a variable are later used to Fourier-transform its corresponding ODE.

```@docs
DifferentialEquation
add_harmonic!
get_variables(::DifferentialEquation)
get_independent_variables(::DifferentialEquation)
```