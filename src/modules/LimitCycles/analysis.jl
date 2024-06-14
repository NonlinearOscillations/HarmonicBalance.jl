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
