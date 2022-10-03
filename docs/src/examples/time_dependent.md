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
@variables ω0, γ, λ, F, x, θ, η, α, ω, t, x(t)

eq =  d(d(x,t),t) + γ*d(x,t) + ω0^2*(1-λ*cos(2*ω*t))*x + α*x^3 + η*d(x,t)*x^2 ~ F*cos(ω*t+θ)

diff_eq = DifferentialEquation(eq, x)
add_harmonic!(diff_eq, x, ω); # single-frequency ansatz

harmonic_eq = get_harmonic_equations(diff_eq);
```
```
A set of 2 harmonic equations
Variables: u1(T), v1(T)
...
```

The object `harmonic_eq` encodes Eq. \eqref{eq:harmeq}.   

We now wish to parse this input into [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) and use its powerful ODE solvers. The desired object here is `DifferentialEquations.ODEProblem`, which is then fed into `DifferentialEquations.solve`.

## Evolving from an initial condition

Given $\mathbf{u}(T_0)$, what is $\mathbf{u}(T)$ at future times?

For constant parameters, a [`HarmonicEquation`](@ref HarmonicBalance.HarmonicEquation) object can be fed into the constructor of [`ODEProblem`](@ref HarmonicBalance.TimeEvolution.ODEProblem). The syntax is similar to DifferentialEquations.jl : 
```julia
import HarmonicBalance.TimeEvolution: ODEProblem, DifferentialEquations.solve
x0 = [0.0; 0.] # initial condition
fixed = (ω0 => 1.0,γ => 1E-2, λ => 5E-2, F => 1E-3,  α => 1., η=>0.3, θ => 0, ω=>1.) # parameter values

ode_problem = ODEProblem(harmonic_eq, fixed, x0 = x0, timespan = (0,1000))
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
plot(time_evo, ["u1", "v1"], harmonic_eq)
```

Running the above code with `x0 = [0., 0.]` and `x0 = [0.2, 0.2]` gives the plots
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/time_dependent/evo_to_steady.png" alignment="center" \>
``` ⠀

Let us compare this to the steady state diagram.
```julia
varied = ω => LinRange(0.9, 1.1, 100)
result = get_steady_states(harmonic_eq, varied, fixed)
plot(result, "sqrt(u1^2 + v1^2)")
```
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/time_dependent/steady.png" alignment="center" \>
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
ode_problem = ODEProblem(harmonic_eq, fixed, sweep=sweep, x0=[0.1;0.0], timespan=(0, 2E4))
time_evo = solve(ode_problem, saveat=100)
plot(time_evo, "sqrt(u1^2 + v1^2)", harmonic_eq)
```
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/time_dependent/sweep_omega.png" alignment="center" \>
``` ⠀
We see the system first evolves from the initial condition towards the low-amplitude steady state. The amplitude increases as the sweep proceeds, with a jump occurring around $\omega = 1.08$ (i.e., time 18000).

Successive sweeps can be combined,
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


## Limit cycles

So far, we have largely focused on finding and analysing steady states, i.e., fixed points of the harmonic equations, which satisfy
```math
\begin{equation} \label{eq:harmeqfull}
\frac{d\mathbf{u}(T)}{dT}  = \bar{\mathbf{F}} (\mathbf{u}) = 0\,.
\end{equation}
```

Fixed points are however merely a subset of possible solutions of Eq. \eqref{eq:harmeqfull} -- strictly speaking, solutions where $\mathbf{u}(T)$ remains time-dependent are allowed. These are quite unusual, since $\bar{\mathbf{F}} (\mathbf{u})$ [is by construction time-independent](@ref intro_hb) and Eq. \eqref{eq:harmeqfull} thus possesses _continuous time-translation symmetry_. The appearance of explicitly time-dependent solutions then consitutes spontaneous time-translation symmetry breaking.

Such solutions, known as _limit cycles_, typically appear as closed periodic trajectories of the harmonic variables $\mathbf{u}(T)$. The simplest way to numerically characterise them is a time-dependent simulation, using a steady-state diagram as a guide.

Here we reconstruct the results of [Zambon et al., Phys Rev. A 102, 023526 (2020)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.102.023526), where limit cycles are shown to appear in a system of two coupled nonlinear oscillators. In this problem, two oscillators $x_1$ and $x_2$, have (the same) damping and Kerr nonlinearity and are linearly coupled,

```math
\begin{align}
\ddot{x}_1+ \gamma \dot{x}_1 + \omega_0^2 x_1 + \alpha x_1^3 + 2J(x_1-x_2) &= F_0 \cos(\omega t) \\
\ddot{x}_2+ \gamma \dot{x}_2 + \omega_0^2 x_2 + \alpha x_2^3 + 2J(x_2-x_1) &= \eta F_0 \cos(\omega t)
\end{align}
```
```julia
using HarmonicBalance
@variables γ, F, α, ω0, F0, η, ω, J, t, x(t), y(t);

diff_eq = DifferentialEquation([d(x,t,2) + γ * d(x,t) + ω0^2 * x + α*x^3+ 2*J*ω0*(x-y) - F0*cos(ω*t), 
            d(y,t,2) + γ * d(y,t) + ω0^2 * y + α*y^3 + 2*J*ω0*(y-x) - η*F0*cos(ω*t)], [x,y])

```

The analysis of Zambon et al. uses a frame rotating at the pump frequency $\omega$ to describe both oscillators. For us, this means we expand both modes using $\omega$ to obtain the harmonic equations.
```julia
add_harmonic!(diff_eq, x, ω)
add_harmonic!(diff_eq, y, ω)

harmonic_eq = get_harmonic_equations(diff_eq)
```
Solving for a range of drive amplitudes $F_0$,
```julia
fixed = (
    ω0 => 1.4504859, # natural frequency of separate modes (in paper's notation, ħω0 - J)
    γ => 27.4E-6,    # damping
    J => 154.1E-6,   # coupling term
    α => 3.867E-7,   # Kerr nonlinearity
    ω => 1.4507941,  # pump frequency, resonant with antisymmetric mode (in paper, ħω0 + J)
    η => -0.08,      # pumping leaking to site 2  (F2 = ηF1)
    F0 => 0.002       # pump amplitude (overriden in sweeps)
)
varied = F0 => LinRange(0.002, 0.03, 50)

result = get_steady_states(harmonic_eq, varied, fixed)
```
```
A steady state result for 50 parameter points

Solution branches:   9
   of which real:    3
   of which stable:  2
```

Let us first see the steady states.
```julia
p1 = plot(result, "u1^2 + v1^2", legend=false)
p2 = plot(result, "u2^2 + v2^2")
plot(p1, p2)
```
```@raw html
<img style="display: block; margin: 0 auto; padding-bottom: 20px" src="../../assets/time_dependent/lc_steady.png" alignment="center"\>
```

According to Zambon et al., a limit cycle solution exists around $F_0 \cong 0.011$, which can be accessed by a jump from branch 1 in an upwards sweep of $F_0$. Since a limit cycle is not a steady state of our harmonic equations, it does not appear in the diagram. We do however see that branch 1 ceases to be stable around $F_0 \cong 0.010$, meaning a jump should occur. 

Let us try and simulate the limit cycle. We could in principle run a time-dependent simulation with a fixed value of $F_0$, but this would require a suitable initial condition. Instead, we will sweep $F_0$ upwards from a low starting value. To observe the dynamics just after the jump has occurred, we follow the sweep by a time interval where the system evolves under fixed parameters.
```julia
import HarmonicBalance.TimeEvolution: ODEProblem, DifferentialEquations.solve
initial_state = result[1][1]

T = 2E6
sweep = ParameterSweep(F0 => (0.002, 0.011), (0,T))

# start from initial_state, use sweep, total time is 2*T
time_problem = ODEProblem(harmonic_eq, initial_state, sweep=sweep, timespan=(0,2*T))
time_evo = solve(time_problem, saveat=100);
```
Inspecting the amplitude as a function of time,
```julia
plot(time_evo, "sqrt(u1^2 + v1^2)", harmonic_eq)
```
```@raw html
<img style="display: block; margin: 0 auto; padding-bottom: 20px" src="../../assets/time_dependent/lc_sweep.png" alignment="center" \>
```

we see that initially the sweep is adiabatic as it proceeds along the steady-state branch 1. At around $T = 2E6$, an instability occurs and $u_1(T)$ starts to rapidly oscillate. At that point, the sweep is stopped. Under free time evolution, the system then settles into a limit-cycle solution where the coordinates move along closed trajectories.

By plotting the $u$ and $v$ variables against each other, we observe the limit cycle shapes in phase space,
```julia
p1 = plot(time_evo, ["u1", "v1"], harmonic_eq)
p2 = plot(time_evo, ["u2", "v2"], harmonic_eq)
plot(p1, p2)
```
```@raw html
<img style="display: block; margin: 0 auto;" src="../../assets/time_dependent/lc_uv.png" alignment="center" \>
```
