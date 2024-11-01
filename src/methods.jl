abstract type HarmonicBalanceMethod end

Base.@kwdef struct TotalDegree <: HarmonicBalanceMethod
    gamma = cis(2π * rand())
    thread::Bool = Threads.nthreads() > 1
    tracker_options = HomotopyContinuation.TrackerOptions()
    endgame_options = HomotopyContinuation.EndgameOptions()
    compile = HomotopyContinuation.COMPILE_DEFAULT[]
    seed = 0xd8e5d8df
end

Base.@kwdef struct Polyhedral <: HarmonicBalanceMethod
    only_non_zero = false
    thread::Bool = Threads.nthreads() > 1
    tracker_options = HomotopyContinuation.TrackerOptions()
    endgame_options = HomotopyContinuation.EndgameOptions()
    compile = HomotopyContinuation.COMPILE_DEFAULT[]
    seed = 0xd8e5d8df
end

Base.@kwdef struct WarmUp <: HarmonicBalanceMethod
    perturbation_size::ComplexF64 = 1e-6 + 1e-6 * im
    index::Union{Int,EndpointRanges.Endpoint} = EndpointRanges.iend ÷ 2
    thread::Bool = Threads.nthreads() > 1
    tracker_options = HomotopyContinuation.TrackerOptions()
    endgame_options = HomotopyContinuation.EndgameOptions()
    compile = HomotopyContinuation.COMPILE_DEFAULT[]
    seed = 0xd8e5d8df
end

thread(method::HarmonicBalanceMethod) = method.thread
tracker(method::HarmonicBalanceMethod) = method.tracker_options
endgame(method::HarmonicBalanceMethod) = method.endgame_options
compile(method::HarmonicBalanceMethod) = method.compile
seed(method::HarmonicBalanceMethod) = method.seed
function alg_default_options(method::HarmonicBalanceMethod)
    return (
        threading=thread(method),
        tracker_options=tracker(method),
        endgame_options=endgame(method),
        compile=compile(method),
        seed=seed(method),
    )
end

alg_specific_options(method::TotalDegree) = (gamma=method.gamma,)
alg_specific_options(method::Polyhedral) = (only_non_zero=method.only_non_zero,)
function alg_specific_options(method::WarmUp)
    return (perturbation_size=method.perturbation_size, index=method.index)
end

function Base.show(io::IO, m::WarmUp)
    println(io, "Warm up method:")
    println(io, "perturbation_size: ", m.perturbation_size)
    println(io, "Threading:         ", thread(m))
    println(io, "Compile:           ", compile(m))
    return println(io, "Seed:              ", seed(m))
end
function Base.show(io::IO, m::TotalDegree)
    println(io, "Total degree method:")
    println(io, "Gamma:     ", m.gamma)
    println(io, "Threading: ", thread(m))
    println(io, "Compile:   ", compile(m))
    return println(io, "Seed:      ", seed(m))
end
function Base.show(io::IO, m::Polyhedral)
    println(io, "Polyhedral method:")
    println(io, "Zero solutions: ", !m.only_non_zero)
    println(io, "Threading:      ", thread(m))
    println(io, "Compile:        ", compile(m))
    return println(io, "Seed:           ", seed(m))
end
