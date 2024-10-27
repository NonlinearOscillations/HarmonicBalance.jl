abstract type HarmonicBalanceMethod end

Base.@kwdef struct TotalDegree <: HarmonicBalanceMethod
    thread::Bool = false
    tracker_options = TrackerOptions()
    endgame_options = EndgameOptions()
    compile = HomotopyContinuation.COMPILE_DEFAULT[]
end

Base.@kwdef struct WarmUp <: HarmonicBalanceMethod
    perturb::ComplexF64 = 1e-6 + 1e-6*im
    index::Union{Int,Symbol} = :middle
    thread::Bool = false
    tracker_options = TrackerOptions()
    endgame_options = EndgameOptions()
    compile = HomotopyContinuation.COMPILE_DEFAULT[]
end
