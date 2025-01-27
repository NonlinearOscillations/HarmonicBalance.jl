function classify_unique!(res::Result, Δω; class_name="unique_cycle")

    # 1st degeneracy: arbitrary sign of Δω
    i1 = _symidx(Δω, res)
    c1 = classify_solutions(res, soln -> _is_physical(soln) && real(soln[i1]) >= 0)

    # 2nd degeneracy: ambiguity in the fixed phase, manifests as the sign of var
    var = HarmonicBalance._remove_brackets(get_cycle_variables(res.problem.eom, Δω)[1])
    i2 = _symidx(var, res)
    c2 = classify_solutions(res, soln -> _is_physical(soln) && real(soln[i2]) >= 0)

    return res.classes[class_name] = map(.*, c1, c2)
end

# if abs(ω_lc) < tol, set all classifications to false
# TOLERANCE HARDCODED FOR NOW
function _classify_limit_cycles!(res::Result, ω_lc::Num)
    ω_lc_idx = findfirst(x -> isequal(x, ω_lc), res.problem.variables)
    for idx in CartesianIndices(res.solutions),
        c in filter(x -> x != "binary_labels", keys(res.classes))

        res.classes[c][idx] .*= abs.(getindex.(res.solutions[idx], ω_lc_idx)) .> 1e-10
    end

    classify_unique!(res, ω_lc)

    unique_stable = find_branch_order(
        map(.*, res.classes["stable"], res.classes["unique_cycle"])
    )

    # branches which are unique but never stable
    unique_unstable = setdiff(
        find_branch_order(map(.*, res.classes["unique_cycle"], res.classes["physical"])),
        unique_stable,
    )
    return order_branches!(res, vcat(unique_stable, unique_unstable)) # shuffle to have relevant branches first
end
