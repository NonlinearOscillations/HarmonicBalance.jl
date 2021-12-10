export add_harmonic!, add_Hopf!
import Symbolics.get_variables


"""
$(TYPEDSIGNATURES)
Add the harmonic `ω` to the harmonic ansatz used to expand the variable `var` in `diff_eom`.
"""
function add_harmonic!(diff_eom::DifferentialEquation, var::Num, ω)
    diff_eom.harmonics[var] = unique(cat(diff_eom.harmonics[var], ω, dims=1))
    diff_eom
end


"""
$(TYPEDSIGNATURES)
Return the dependent variables of `diff_eom`.
"""
get_variables(diff_eom::DifferentialEquation) = collect(keys(diff_eom.equations))


is_harmonic(diff_eom::DifferentialEquation, t::Num) = all([is_harmonic(eq, t) for eq in values(diff_eom.equations)])

"Pretty printing of the newly defined types" 
function show_fields(object)
    for field in fieldnames(typeof(object)) # display every field
        display(string(field)); display(getfield(object, field))
    end
end

"""
$(TYPEDSIGNATURES)
Return the independent dependent variables of `diff_eom`.
"""
function get_independent_variables(diff_eom::DifferentialEquation)
    Num.(flatten(unique([x.val.arguments for x in keys(diff_eom.equations)])))
end

show(eom::DifferentialEquation) = show_fields(eom)


