import HarmonicBalance.load

# load the previously-saved result
current_path = @__DIR__
@test load(current_path * "/parametron_result.jld2") isa HarmonicBalance.Result skip = true
