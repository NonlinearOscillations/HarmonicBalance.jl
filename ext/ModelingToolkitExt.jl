module ModelingToolkitExt

export ODESystem, ODEProblem

using HarmonicBalance: HarmonicEquation, is_rearranged,
                       rearrange_standard, get_variables, simplify, ParameterList
using ModelingToolkit
import ModelingToolkit: ODESystem, ODEProblem, NonlinearProblem, SteadyStateProblem

swapsides(eq::Equation) = Equation(eq.rhs, eq.lhs)

function declare_parameter(var::Num)
    var_sym = Symbol(var)
    new_var = @parameters $var_sym
    @eval($(var_sym)=first($new_var)) # store the variable under "name" in this namespace
    return eval(var_sym)
end

function ODESystem(eom::HarmonicEquation)
    if !is_rearranged(eom) # check if time-derivatives of the variable are on the right hand side
        eom = rearrange_standard(eom)
    end

    slow_time = (@variables T; T)
    par_names = declare_parameter.(eom.parameters)
    vars = get_variables(eom)

    eqs = deepcopy(eom.equations)
    eqs = swapsides.(eqs)
    eqs = simplify.(expand.(eqs))

    # compute jacobian for performance
    @named sys = ODESystem(eqs, slow_time, vars, par_names)
    return sys
end

function ODEProblem(eom::HarmonicEquation, u0, tspan::Tuple,
        p::ParameterList; in_place = true, kwargs...)
    sys = ODESystem(eom)
    param = ModelingToolkit.varmap_to_vars(p, parameters(sys))
    if !in_place # out-of-place
        prob = ODEProblem{false}(sys, u0, tspan, param; jac = true, kwargs...)
    else # in-place
        prob = ODEProblem{true}(sys, u0, tspan, param; jac = true, kwargs...)
    end
    return prob
end

function NonlinearProblem(
        eom::HarmonicEquation, u0, p::ParameterList; in_place = true, kwargs...)
    ss_prob = SteadyStateProblem(eom, u0, p::ParameterList; in_place = in_place, kwargs...)
    return NonlinearProblem(ss_prob) # gives warning of some internal deprication
end

function SteadyStateProblem(
        eom::HarmonicEquation, u0, p::ParameterList; in_place = true, kwargs...)
    sys = ODESystem(eom)
    param = ModelingToolkit.varmap_to_vars(p, parameters(sys))
    if !in_place # out-of-place
        prob = SteadyStateProblem{false}(sys, u0, param; jac = true, kwargs...)
    else # in-place
        prob = SteadyStateProblem{true}(sys, u0, param; jac = true, kwargs...)
    end
    return prob
end

end # module
