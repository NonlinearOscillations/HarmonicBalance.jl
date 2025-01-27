const StateDict = OrderedDict;
const Solutions(T) = VecOrMat{Vector{Vector{T}}};
const SteadyState(T) = Vector{T};

solution_type(sol::Solutions(T)) where {T} = T

const JacobianFunction(T) = FunctionWrapper{Matrix{T},Tuple{Vector{T}}}

"""

Represents a sweep of one or more parameters of a `HarmonicEquation`.
During a sweep, the selected parameters vary linearly over some timespan and are constant elsewhere.

Sweeps of different variables can be combined using `+`.

# Fields
$(TYPEDFIELDS)

## Examples
```julia-repl
# create a sweep of parameter a from 0 to 1 over time 0 -> 100
julia> @variables a,b;
julia> sweep = AdiabaticSweep(a => [0., 1.], (0, 100));
julia> sweep[a](50)
0.5
julia> sweep[a](200)
1.0

# do the same, varying two parameters simultaneously
julia> sweep = AdiabaticSweep([a => [0.,1.], b => [0., 1.]], (0,100))
```
Successive sweeps can be combined,
```julia-repl
sweep1 = AdiabaticSweep(ω => [0.95, 1.0], (0, 2e4))
sweep2 = AdiabaticSweep(λ => [0.05, 0.01], (2e4, 4e4))
sweep = sweep1 + sweep2
```
multiple parameters can be swept simultaneously,
```julia-repl
sweep = AdiabaticSweep([ω => [0.95;1.0], λ => [5e-2;1e-2]], (0, 2e4))
```
and custom sweep functions may be used.
```julia-repl
ωfunc(t) = cos(t)
sweep = AdiabaticSweep(ω => ωfunc)
```

"""
struct AdiabaticSweep
    """Maps each swept parameter to a function."""
    functions::Dict{Num,Function}

    AdiabaticSweep(functions...) = new(Dict(functions...))
    AdiabaticSweep() = AdiabaticSweep([])
end
