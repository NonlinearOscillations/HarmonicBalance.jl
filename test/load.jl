import HarmonicBalance.load
# using JLD2
import JLD2.JLDFile, JLD2.RelOffset, JLD2.NULL_REFERENCE
import JLD2.track_weakref_if_untracked!

function track_weakref_if_untracked!(f::JLDFile, header_offset::RelOffset, @nospecialize v)
    if header_offset !== NULL_REFERENCE
        if !haskey(f.jloffset, header_offset) || isnothing(f.jloffset[header_offset].value)
            f.jloffset[header_offset] = WeakRef(v)
        end
    end
    nothing
end

# load the previously-saved result
current_path = @__DIR__
@test load(current_path * "/parametron_result.jld2") isa HarmonicBalance.Result
