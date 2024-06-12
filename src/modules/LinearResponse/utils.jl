"Return the indices of u,v pairs belonging to the same natural variable and harmonic from hvars."
function _get_uv_pairs(hvars::Vector{HarmonicVariable})
    u_idx = findall(x -> x.type == "u", hvars)
    pairs = []
    for i in u_idx
        j = findall(
            x ->
                isequal(x.ω, hvars[i].ω) &&
                    x.type == "v" &&
                    isequal(hvars[i].natural_variable, x.natural_variable),
            hvars,
        )
        j = length(j) != 1 && error("no v coordinate found for ", hvars[i]) || j[1]
        push!(pairs, [i, j])
    end
    return pairs
end
