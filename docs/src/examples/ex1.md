# Example 1

some text

markdown: 
```julia
using HarmonicBalance.jl
```

First, we need to define the symbolic variables to be used. The syntax here is inherited from Symbolics.jl
```julia
@variables
```

These are used to define an equation of motion
```julia
dEOM = HarmonicBalance.DifferentialEquation
```

This object now holds the equations and identifies its variables.

Now we must specify which harmonics will be used to expand which variable, let's start with using the single
```julia
HarmonicBalance.add_harmonic!(dEOM, x, ω)
```

The problem is now fully specified - we may convert the differential equation into a 
```julia
averagedEOM
```

We may now inspect the generated equations. The conversion has generated new variables describing each harmonic.
```julia
averagedEOM.equations
```

This object is the central point of the HarmonicBalance.jl framework. It may be fed into an algebraic or a time-dependent solver. 

```julia
averagedEOM
```

We are now ready to use the homotopy continuation method to obtain all the solutions of our algebraic equations. This is typically done across a set of parameters.

```julia
ParameterRange
get_steady_states
```


So far, the equation have been relatively uncomplicated. However, due to the presence of the cubic nonlinearity, the ansatz ...      $`3ω`$ does not fully describe the behaviour of our system. Frequency conversion - to first order, from ω to 3ω - will occur. To see the effect of this, we may expand the motion in both harmonics.
```julia
HarmonicBalance.add_harmonic!(dEOM, x, [ω, 3ω])
```

The generated equations are not perturbative - all harmonics are treated on the same footing. The newly generated variables are (pretty print of VDP types). Proceeding to a solution diagram as before, we may now plot the two harmonics separately.


We see that, although the third harmonic is clearly present, our solution landscape did not change much. In later examples, we will see this is not always the case. 



To sort: hyperlinks
