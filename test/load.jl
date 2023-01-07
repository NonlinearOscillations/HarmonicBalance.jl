import HarmonicBalance.load

# load the previously-saved result
parametron_result_path = joinpath(@__DIR__, "parametron_result.jld2")
@test load(parametron_result_path) isa HarmonicBalance.Result
