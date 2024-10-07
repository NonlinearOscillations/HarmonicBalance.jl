"""
classify_unique!(res::Result, Δω; class_name="unique_cycle")

Classifies solutions in `res` based on the degeneracies related to the sign of `Δω` and the ambiguity in the fixed phase, returning a unique classification.

### Description
This function removes two types of degeneracies in the solutions:
1. **First degeneracy**: The arbitrary sign of `Δω`, classified by checking if the real part of `Δω` is non-negative.
2. **Second degeneracy**: Ambiguity in the fixed phase, which manifests as the sign of a specific variable in the system's equation of motion (EOM). This variable is extracted from the cycle of the system and is also classified based on its real part being non-negative.

The combined classification, a logical AND (`.*`) of the two degeneracies, is stored in the `res.classes` dictionary under the key `class_name`.
"""
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
