import HarmonicBalance.classify_solutions

function classify_unique!(res::Result, Δω; class_name="unique")

    # 1st degeneracy: arbitrary sign of Δω
    c1 = classify_solutions(res, string(Δω) * ">= 0")

    # 2nd degeneracy: ambiguity in the fixed phase, manifests as the sign of var
    var = HarmonicBalance._remove_brackets(get_cycle_variables(res.problem.eom, Δω)[1])
    c2 = classify_solutions(res, string(var) * ">= 0") 

    res.classes[class_name] = map(.*, c1, c2)
end