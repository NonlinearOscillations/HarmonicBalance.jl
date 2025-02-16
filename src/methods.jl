"""
    HarmonicBalanceMethod

Abstract type for harmonic balance methods.
"""
abstract type HarmonicBalanceMethod end

"""
    TotalDegree

The Total Degree homotopy method performs a homotopy ``H(x, t) = γ t G(x) + (1-t) F(x)``
from the trivial polynomial system ``F(x) =xᵢ^{dᵢ} +aᵢ`` with the maximal degree ``dᵢ``
determined by the [Bezout bound](https://en.wikipedia.org/wiki/B%C3%A9zout%27s_theorem).
The method guarantees to find all solutions, however, it comes with a high computational cost.
See [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/guides/totaldegree/)
for more information.

# Fields
$(TYPEDFIELDS)
"""
struct TotalDegree{T<:Complex} <: HarmonicBalanceMethod
    """Complex multiplying factor of the start system G(x) for the homotopy"""
    gamma::T
    """Boolean indicating if threading is enabled."""
    thread::Bool
    """Options for the tracker."""
    tracker_options::HomotopyContinuation.TrackerOptions
    """Options for the endgame."""
    endgame_options::HomotopyContinuation.EndgameOptions
    """Compilation options."""
    compile::Union{Bool,Symbol}
    """Seed for random number generation."""
    seed::UInt32

    function TotalDegree(;
        gamma::T=cis(2π * rand(Random.MersenneTwister(0xd8e5d8df))),
        thread::Bool=Threads.nthreads() > 1,
        tracker_options::HomotopyContinuation.TrackerOptions=HomotopyContinuation.TrackerOptions(),
        endgame_options::HomotopyContinuation.EndgameOptions=HomotopyContinuation.EndgameOptions(),
        compile::Union{Bool,Symbol}=HomotopyContinuation.COMPILE_DEFAULT[],
        seed::UInt32=0xd8e5d8df,
    ) where {T<:Complex}
        return new{T}(gamma, thread, tracker_options, endgame_options, compile, seed)
    end
    function TotalDegree{T}(;
        gamma::Complex=cis(2π * rand(Random.MersenneTwister(0xd8e5d8df), T)),
        thread::Bool=Threads.nthreads() > 1,
        tracker_options::HomotopyContinuation.TrackerOptions=HomotopyContinuation.TrackerOptions(),
        endgame_options::HomotopyContinuation.EndgameOptions=HomotopyContinuation.EndgameOptions(),
        compile::Union{Bool,Symbol}=HomotopyContinuation.COMPILE_DEFAULT[],
        seed::UInt32=0xd8e5d8df,
    ) where {T<:Complex}
        return new{T}(gamma, thread, tracker_options, endgame_options, compile, seed)
    end
end

"""
    Polyhedral

The Polyhedral homotopy method constructs a homotopy based on the polyhedral
structure of the polynomial system. It is more efficient than the Total Degree method for
sparse systems, meaning most of the coefficients are zero. It can be especially useful if
you don't need to find the zero solutions (`only_non_zero = true`), resulting in a speed up.
See [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/guides/polyhedral/)
for more information.

# Fields
$(TYPEDFIELDS)
"""
struct Polyhedral{T} <: HarmonicBalanceMethod
    """Boolean indicating if only non-zero solutions are considered."""
    only_non_zero::Bool
    """Boolean indicating if threading is enabled."""
    thread::Bool
    """Options for the tracker."""
    tracker_options::HomotopyContinuation.TrackerOptions
    """Options for the endgame."""
    endgame_options::HomotopyContinuation.EndgameOptions
    """Compilation options."""
    compile::Union{Bool,Symbol}
    """Seed for random number generation."""
    seed::UInt32

    function Polyhedral(;
        only_non_zero::Bool=false,
        thread::Bool=Threads.nthreads() > 1,
        tracker_options::HomotopyContinuation.TrackerOptions=HomotopyContinuation.TrackerOptions(),
        endgame_options::HomotopyContinuation.EndgameOptions=HomotopyContinuation.EndgameOptions(),
        compile::Union{Bool,Symbol}=HomotopyContinuation.COMPILE_DEFAULT[],
        seed::UInt32=0xd8e5d8df,
    )
        return new{ComplexF64}(
            only_non_zero, thread, tracker_options, endgame_options, compile, seed
        )
    end
    function Polyhedral{T}(;
        only_non_zero::Bool=false,
        thread::Bool=Threads.nthreads() > 1,
        tracker_options::HomotopyContinuation.TrackerOptions=HomotopyContinuation.TrackerOptions(),
        endgame_options::HomotopyContinuation.EndgameOptions=HomotopyContinuation.EndgameOptions(),
        compile::Union{Bool,Symbol}=HomotopyContinuation.COMPILE_DEFAULT[],
        seed::UInt32=0xd8e5d8df,
    ) where {T<:Complex}
        return new{T}(
            only_non_zero, thread, tracker_options, endgame_options, compile, seed
        )
    end
end

"""
    WarmUp

The Warm Up method prepares a warmup system with the Total Degree method using the parameter
at `index` perturbed by `perturbation_size`. The warmup system is used to perform a homotopy
using all other systems in the parameter sweep. It is very efficient for systems with minimal
bifurcation in the parameter sweep. The Warm Up method should in theory guarantee to find all
solutions, however, if the `start_parameters` is not proper (to close to the real line)
it could miss some solutions.

See[HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/guides/many-systems/)
for more information.

# Fields
$(TYPEDFIELDS)
"""
struct WarmUp{T} <: HarmonicBalanceMethod
    """Method used for the warmup system."""
    warm_up_method::Union{TotalDegree{T},Polyhedral{T}}
    """Start parameters."""
    start_parameters::Vector{T}
    """Boolean indicating if threading is enabled."""
    thread::Bool
    """Check if zero is a root"""
    check_zero::Bool
    """Options for the tracker."""
    tracker_options::HomotopyContinuation.TrackerOptions
    """Options for the endgame."""
    endgame_options::HomotopyContinuation.EndgameOptions
    """Compilation options."""
    compile::Union{Bool,Symbol}
    """Seed for random number generation."""
    seed::UInt32

    function WarmUp(;
        warm_up_method::Union{TotalDegree,Polyhedral}=TotalDegree(),
        start_parameters::Vector=Vector{ComplexF64}(),
        thread::Bool=Threads.nthreads() > 1,
        check_zero::Bool=true,
        tracker_options::HomotopyContinuation.TrackerOptions=HomotopyContinuation.TrackerOptions(),
        endgame_options::HomotopyContinuation.EndgameOptions=HomotopyContinuation.EndgameOptions(),
        compile::Union{Bool,Symbol}=HomotopyContinuation.COMPILE_DEFAULT[],
        seed::UInt32=0xd8e5d8df,
    )
        T1 = eltype(start_parameters)
        T2 = typeof(warm_up_method).parameters[1]
        if T1 != T2
            @warn "The start parameters and the method parameters do not match.
            The start parameters are converted to the method parameters."
            start_parameters = Vector{T2}(start_parameters)
        end
        return new{T2}(
            warm_up_method,
            start_parameters,
            thread,
            check_zero,
            tracker_options,
            endgame_options,
            compile,
            seed,
        )
    end
    function WarmUp{T}(;
        warm_up_method::Union{TotalDegree,Polyhedral}=Polyhedral{T}(),
        start_parameters::Vector=Vector{T}(),
        thread::Bool=Threads.nthreads() > 1,
        check_zero::Bool=true,
        tracker_options::HomotopyContinuation.TrackerOptions=HomotopyContinuation.TrackerOptions(),
        endgame_options::HomotopyContinuation.EndgameOptions=HomotopyContinuation.EndgameOptions(),
        compile::Union{Bool,Symbol}=HomotopyContinuation.COMPILE_DEFAULT[],
        seed::UInt32=0xd8e5d8df,
    ) where {T<:Complex}
        return new{T}(
            warm_up_method,
            start_parameters,
            thread,
            check_zero,
            tracker_options,
            endgame_options,
            compile,
            seed,
        )
    end
end

"""
    thread(method::HarmonicBalanceMethod) -> Bool

Returns whether threading is enabled for the given method. The number of available threads
is controlled by the environment variable `JULIA_NUM_THREADS`.
"""
thread(method::HarmonicBalanceMethod) = method.thread

"""
    tracker(method::HarmonicBalanceMethod) -> HomotopyContinuation.TrackerOptions

Returns the tracker options for the given method. See `HomotopyContinuation.TrackerOptions`
for the available options.
"""
tracker(method::HarmonicBalanceMethod) = method.tracker_options

"""
    endgame(method::HarmonicBalanceMethod) -> HomotopyContinuation.EndgameOptions

Returns the endgame options for the given method. See `HomotopyContinuation.EndgameOptions`
for the available options.
"""
endgame(method::HarmonicBalanceMethod) = method.endgame_options

"""
    compile(method::HarmonicBalanceMethod) -> Union{Bool,Symbol}

Returns the compile options for the given method. If `true` then a system is compiled to a
straight line program for evaluation. This induces a compilation overhead. If `false` then
the generated program is only interpreted. This is slower than the compiled version,
but does not introduce compilation overhead.
"""
compile(method::HarmonicBalanceMethod) = method.compile

"""
    seed(method::HarmonicBalanceMethod) -> UInt32

Returns the seed for random number generation for the given method.
"""
seed(method::HarmonicBalanceMethod) = method.seed

"""
    alg_default_options(method::HarmonicBalanceMethod) -> NamedTuple

Returns a named tuple of default algorithm options for the given method.
"""
function alg_default_options(method::HarmonicBalanceMethod)
    return (
        threading=thread(method),
        tracker_options=tracker(method),
        endgame_options=endgame(method),
        compile=compile(method),
        seed=seed(method),
    )
end

"""
    alg_specific_options(method::TotalDegree) -> NamedTuple

Returns a named tuple of specific algorithm options for the Total Degree method.
"""
alg_specific_options(method::TotalDegree) = (gamma=method.gamma,)

"""
    alg_specific_options(method::Polyhedral) -> NamedTuple

Returns a named tuple of specific algorithm options for the Polyhedral method.
"""
alg_specific_options(method::Polyhedral) = (only_non_zero=method.only_non_zero,)

"""
    method_symbol(m::Polyhedral) -> Symbol

Returns the symbol for the Polyhedral method identified by HomotopyContinuation/jl.
"""
method_symbol(m::Polyhedral) = :polyhedral

"""
    method_symbol(m::TotalDegree) -> Symbol

Returns the symbol for the Total Degree method identified by HomotopyContinuation.jl.
"""
method_symbol(m::TotalDegree) = :total_degree

"""
    Base.show(io::IO, m::WarmUp)

Displays information about the Warm Up method.
"""
function Base.show(io::IO, m::WarmUp)
    println(io, "$(typeof(m)) method:")
    println(io, "Threading:         ", thread(m))
    println(io, "Compile:           ", compile(m))
    return println(io, "Seed:              ", seed(m))
end

"""
    Base.show(io::IO, m::TotalDegree)

Displays information about the Total Degree method.
"""
function Base.show(io::IO, m::TotalDegree)
    println(io, "$(typeof(m)) method:")
    println(io, "Gamma:     ", m.gamma)
    println(io, "Threading: ", thread(m))
    println(io, "Compile:   ", compile(m))
    return println(io, "Seed:      ", seed(m))
end

"""
    Base.show(io::IO, m::Polyhedral)

Displays information about the Polyhedral method.
"""
function Base.show(io::IO, m::Polyhedral)
    println(io, "$(typeof(m)) method:")
    println(io, "Zero solutions: ", !m.only_non_zero)
    println(io, "Threading:      ", thread(m))
    println(io, "Compile:        ", compile(m))
    return println(io, "Seed:           ", seed(m))
end
