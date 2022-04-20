# Introduction: time-dependent simulations

Most of HarmonicBalance.jl is focused on finding and analysing the steady states. Such states contain no information about transient behaviour, which is crucial to answer the following.

* Given an initial condition, which steady state does the system evolve into?
* How does the system behave if its parameters are varied in time?

It is straightforward to evolve the full equation of motion using an ODE solver. However, tracking oscillatory behaviour is computationally expensive.

In the [background](@ref intro_hb), we showed that nonlinear driven systems may be reduced to harmonic equations

```math
\begin{equation} \label{eq:harmeq}
\frac{d\mathbf{u}(T)}{dT}  = \bar{\mathbf{F}} (\mathbf{u})\,,
\end{equation}
```

As long as the chosen harmonics constituting $\mathbf{u}(T)$ capture the system's behaviour, we may numerically evolve  Eq. \eqref{eq:harmeq} instead of the full problem. Since the components of $\mathbf{u}(T)$ only vary very slowly (and are constant in a steady state), this is usually _vastly_ more efficient than evolving the full problem. 

Here we primarily demonstrate on the [parametrically driven oscillator](@ref parametron). The relevant notebooks may be found [in the example repo](https://github.com/NonlinearOscillations/HarmonicBalance-notebooks).

We start by defining our system.
```julia
using HarmonicBalance
@variables Ω, γ, λ, F, x, θ, η, α, ω, ψ, T, t, x(t)

natural_equation =  d(d(x,t),t) + γ*d(x,t) + Ω^2*(1-λ*cos(2*ω*t+ψ))*x + α*x^3 + η*d(x,t)*x^2
force =  F*cos(ω*t+θ)
dEOM = HarmonicBalance.DifferentialEquation(natural_equation + force, x)
HarmonicBalance.add_harmonic!(dEOM, x, ω); # single-frequency ansatz

# construct the harmonic equations
harmonic_eqs = HarmonicBalance.get_harmonic_equations(dEOM)
```
```
A set of 2 harmonic equations
Variables: u1(T), v1(T)
...
```

The object `harmonic_eqs` encodes Eq. \eqref{eq:harmeq}.   

**We now wish to parse this input into [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) and use its powerful ODE solvers.** The desired object here is `DifferentialEquations.ODEProblem`, which is then fed into `DifferentialEquations.solve`.

## Evolving from an initial condition

Given $\mathbf{u}(T_0)$, what is $\mathbf{u}(T)$ at future times?

For constant parameters, a [`HarmonicEquation`](@ref HarmonicBalance.HarmonicEquation) object can be fed into the constructor of [`ODEProblem`](@ref HarmonicBalance.TimeEvolution.ODEProblem). The syntax is similar to DifferentialEquations.jl : 
```julia
import HarmonicBalance.TimeEvolution: ODEProblem, DifferentialEquations.solve
x0 = [0.0; 0.] # initial condition
fixed = (Ω => 1.0,γ => 1E-2, λ => 5E-2, F => 1E-3,  α => 1., η=>0.3, θ => 0, ψ => 0, ω=>1.) # parameter values

ode_problem = ODEProblem(harmonic_eqs, fixed, x0 = x0, timespan = (0,1000))
```

```
ODEProblem with uType Vector{Float64} and tType Int64. In-place: true
timespan: (0, 1000)
u0: 2-element Vector{Float64}:
 0.0
 0.0
```

DifferentialEquations.jl takes it from here - we only need to use `solve`.

```julia
time_evo = solve(ode_problem, saveat=1.);
plot(getindex.(time_evo.u, 1), getindex.(time_evo.u,2))
```

Running the above code with `x0 = [0., 0.]` and `x0 = [0.2, 0.2]` gives the plots
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/time_dependent/ss_approach.png" width="900" alignment="center" \>
``` ⠀

Let us compare this to the steady state diagram.
```julia
range = ω => LinRange(0.9, 1.1, 100)           # range of parameter values
solutions = get_steady_states(harmonic_eqs, range, fixed)
plot_1D_solutions(solutions, x="ω", y="sqrt(u1^2 + v1^2)*sign(u1)");
```
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/time_dependent/parametron_response.png" width="600" alignment="center" \>
``` ⠀

Clearly when evolving from `x0 = [0.,0.]`, the system ends up in the low-amplitude branch 2. With `x0 = [0.2, 0.2]`, the system ends up in branch 3.

## Parameter sweeps

Experimentally, the primary means of exploring the steady state landscape is an adiabatic sweep one or more of the system parameters. This takes the system along a solution branch. If this branch disappears or becomes unstable, a jump occurs.

The object [`ParameterSweep`](@ref HarmonicBalance.TimeEvolution.ParameterSweep) specifies a sweep, which is then used as an optional `sweep` keyword in the `ODEProblem` constructor.
```julia
sweep = ParameterSweep(ω => (0.9,1.1), (0, 2E4))
```
The sweep linearly interpolates between $\omega = 0.9$ at time 0 and $\omega  = 1.1$ at time 2e4. For earlier/later times, $\omega$ is constant.

Let us now define a new `ODEProblem` which incorporates `sweep` and again use `solve`:
```julia
ode_problem = ODEProblem(harmonic_eqs, fixed, sweep=sweep, x0=[0.1;0.0], timespan=(0, 2E4))
time_evo = solve(prob, saveat=100)
plot(time_evo.t, sqrt.(getindex.(time_evo.u,1).^2 .+ getindex.(time_evo,2).^2))
```
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/time_dependent/sweep_omega.png" width="450" alignment="center" \>
``` ⠀
We see the system first evolves from the initial condition towards the low-amplitude steady state. The amplitude increases as the sweep proceeds, with a jump occurring around $\omega = 1.08$ (i.e., time 18000).

Successive weeps can be combined,
```julia
sweep1 = ParameterSweep(ω => [0.95, 1.0], (0, 2E4))
sweep2 = ParameterSweep(λ => [0.05, 0.01], (2E4, 4E4))
sweep = sweep1 + sweep2
```
multiple parameters can be swept simultaneously,
```julia
sweep = ParameterSweep([ω => [0.95;1.0], λ => [5E-2;1E-2]], (0, 2E4))
```

and custom sweep functions may be used.
```julia
ωfunc(t) = cos(t)
sweep = ParameterSweep(ω => ωfunc)
```

Either of the above can be fed into `ODEProblem` and simulated. See [here](https://github.com/NonlinearOscillations/HarmonicBalance-notebooks) for more examples.