
# Krylov-Bogoliubov Averaging Method {#Krylov-Bogoliubov}

The Krylov-Bogoliubov averaging method is an alternative high-frequency expansion technique used to analyze dynamical systems. Unlike the [Harmonic Balance method](https://en.wikipedia.org/wiki/Harmonic_balance), which is detailed in the [background section](/background/harmonic_balance#intro_hb), the Krylov-Bogoliubov method excels in computing higher orders in $1/\omega$, enabling the capture of faster dynamics within a system.

## Purpose and Advantages {#Purpose-and-Advantages}

The primary advantage of the Krylov-Bogoliubov method lies in its ability to delve deeper into high-frequency components, allowing a more comprehensive understanding of fast dynamical behaviors. By leveraging this technique, one can obtain higher-order approximations that shed light on intricate system dynamics.

However, it&#39;s essential to note a limitation: this method cannot handle multiple harmonics within a single variable, unlike some other high-frequency expansion methods.

## Usage

To compute the Krylov-Bogoliubov averaging method within your system, utilize the function `get_krylov_equations`. This function is designed specifically to implement the methodology and derive the equations necessary to analyze the system dynamics using this technique.

### Function Reference {#Function-Reference}
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.KrylovBogoliubov.get_krylov_equations-manual-Krylov-Bogoliubov_method' href='#HarmonicBalance.KrylovBogoliubov.get_krylov_equations-manual-Krylov-Bogoliubov_method'><span class="jlbinding">HarmonicBalance.KrylovBogoliubov.get_krylov_equations</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
get_krylov_equations(
    diff_eom::DifferentialEquation;
    order,
    fast_time,
    slow_time
)

```


Apply the Krylov-Bogoliubov averaging method to a specific `order` to obtain a set of ODEs (the slow-flow equations) governing the harmonics of `diff_eom`.

The harmonics evolve in `slow_time`, the oscillating terms themselves in `fast_time`. If no input is used, a variable T is defined for `slow_time` and `fast_time` is taken as the independent variable of `diff_eom`.

Krylov-Bogoliubov averaging method can be applied up to `order = 2`.

**Example**

```julia
julia> @variables t, x(t), ω0, ω, F;

# enter the simple harmonic oscillator
julia> diff_eom = DifferentialEquation( d(x,t,2) + ω0^2 * x ~ F *cos(ω*t), x);

# expand x in the harmonic ω
julia> add_harmonic!(diff_eom, x, ω);

# get equations for the harmonics evolving in the slow time T to first order
julia> harmonic_eom = get_krylov_equations(diff_eom, order = 1)

A set of 2 harmonic equations
Variables: u1(T), v1(T)
Parameters: ω, F, ω0

Harmonic ansatz:
xˍt(t) =
x(t) = u1(T)*cos(ωt) + v1(T)*sin(ωt)

Harmonic equations:

((1//2)*(ω^2)*v1(T) - (1//2)*(ω0^2)*v1(T)) / ω ~ Differential(T)(u1(T))

((1//2)*(ω0^2)*u1(T) - (1//2)*F - (1//2)*(ω^2)*u1(T)) / ω ~ Differential(T)(v1(T))
```



[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/modules/KrylovBogoliubov.jl#L25)

</details>


For further information and a detailed understanding of this method, refer to [Krylov-Bogoliubov averaging method on Wikipedia](https://en.wikipedia.org/wiki/Krylov%E2%80%93Bogoliubov_averaging_method).
