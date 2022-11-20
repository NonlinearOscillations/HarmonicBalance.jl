using Symbolics
using OrderedCollections
import HomotopyContinuation

export DifferentialEquation, HarmonicVariable, HarmonicEquation,Problem, Result

const ParameterRange = OrderedDict{Num, Vector{Float64}}; export ParameterRange
const ParameterList = OrderedDict{Num, Float64}; export ParameterList;
const StateDict = OrderedDict{Num, ComplexF64}; export StateDict
const SteadyState = Vector{ComplexF64}; export SteadyState;
const ParameterVector = Vector{Float64}; export ParameterVector;

function Base.getindex(p::ParameterRange, idx::Int...)
   lengths = [length(a) for a in values(p)]
   indices = CartesianIndices(Tuple(lengths))[idx...]
   return [val[indices[j]] for (j,val) in enumerate(values(p))]
end

"""
$(TYPEDEF)

Holds differential equation(s) of motion and a set of harmonics to expand each variable.
This is the primary input for `HarmonicBalance.jl` ; after inputting the equations, the harmonics
    ansatz needs to be specified using `add_harmonic!`.

# Fields
$(TYPEDFIELDS)

## Example
```julia-repl
julia> @variables t, x(t), y(t), ω0, ω, F, k;

# equivalent ways to enter the simple harmonic oscillator
julia> DifferentialEquation(d(x,t,2) + ω0^2 * x - F * cos(ω*t), x);
julia> DifferentialEquation(d(x,t,2) + ω0^2 * x ~ F * cos(ω*t), x);

# two coupled oscillators, one of them driven
julia> DifferentialEquation([d(x,t,2) + ω0^2 * x - k*y, d(y,t,2) + ω0^2 * y - k*x] .~ [F * cos(ω*t), 0], [x,y]);
```
"""
mutable struct DifferentialEquation
    """Assigns to each variable an equation of motion."""
    equations::OrderedDict{Num, Equation}
    """Assigns to each variable a set of harmonics."""
    harmonics::OrderedDict{Num, Vector{Num}}

    DifferentialEquation(eqs) = new(eqs, OrderedDict([[var, []] for var in keys(eqs)]))

    # uses the above constructor if no harmonics defined
    DifferentialEquation(eqs::Vector{Equation}, vars::Vector{Num}) = DifferentialEquation(OrderedDict(zip(vars, eqs)))

    # if expressions are entered instead of equations, automatically set them = 0
    DifferentialEquation(exprs::Vector{Num}, vars::Vector{Num}) = DifferentialEquation(OrderedDict(zip(vars, exprs .~ zeros(Int, length(exprs)))))

    DifferentialEquation(arg1, arg2) = DifferentialEquation(arg1 isa Vector ? arg1 : [arg1], arg2 isa Vector ? arg2 : [arg2])
end


function show(io::IO, diff_eq::DifferentialEquation)
    println(io, "System of ", length(keys(diff_eq.equations)), " differential equations")
    println(io, "Variables:       ", join(keys(diff_eq.equations), ", "))
    print(io, "Harmonic ansatz: ")
    for var in keys(diff_eq.harmonics)
        print(io, string(var), " => ", join(string.(diff_eq.harmonics[var]), ", "))
        print(io, ";   ")
    end
    println(io, "\n")
    [println(io, eq) for eq in values(diff_eq.equations)]
end


"""
$(TYPEDEF)

Holds a variable stored under `symbol` describing the harmonic `ω` of `natural_variable`.

# Fields
$(TYPEDFIELDS)
"""
mutable struct HarmonicVariable
    """Symbol of the variable in the HarmonicBalance namespace."""
    symbol::Num
    """Human-readable labels of the variable, used for plotting."""
    name::String
    """Type of the variable (u or v for quadratures, a for a constant, Hopf for Hopf etc.)"""
    type::String
    """The harmonic being described."""
    ω::Num
    """The natural variable whose harmonic is being described."""
    natural_variable::Num
end


function show(io::IO, hv::HarmonicVariable)
    println(io, "Harmonic variable ", string.(hv.symbol) * " for harmonic ", string(hv.ω), " of ", string(hv.natural_variable))
end


"""
$(TYPEDEF)

Holds a set of algebraic equations governing the harmonics of a `DifferentialEquation`.

# Fields
$(TYPEDFIELDS)
"""
mutable struct HarmonicEquation
    """A set of equations governing the harmonics."""
    equations::Vector{Equation}
    """A set of variables describing the harmonics."""
    variables::Vector{HarmonicVariable}
    """The parameters of the equation set."""
    parameters::Vector{Num}
    "The natural equation (before the harmonic ansatz was used)."
    natural_equation::DifferentialEquation

    # use a self-referential constructor with _parameters
    HarmonicEquation(equations, variables, nat_eq) = (x = new(equations, variables, Vector{Num}([]), nat_eq); x.parameters=_parameters(x); x)
end


function show(io::IO, eom::HarmonicEquation)
    println(io, "A set of ", length(eom.equations), " harmonic equations")
    println(io, "Variables: ", join(string.(get_variables(eom)), ", "))
    println(io, "Parameters: ", join(string.(eom.parameters), ", "))
    println(io, "\nHarmonic ansatz: ", _show_ansatz(eom))
    println(io, "\nHarmonic equations:")
    [println(io, "\n", eq) for eq in eom.equations]
end

"""Gives the relation between `var` and the underlying natural variable."""
function _show_ansatz(var::HarmonicVariable)
    t = var.natural_variable.val.arguments
    t = length(t) == 1 ? string(t[1]) : error("more than 1 independent variable")
    ω = string(var.ω)
    terms = Dict("u" => "*cos(" * ω * t * ")", "v" => "*sin(" * ω * t * ")", "a" => "")
    string(string(var.symbol) * terms[var.type])
end


"""Gives the full harmonic ansatz used to construct `eom`."""
function _show_ansatz(eom::HarmonicEquation)
    output = ""
    for nat_var in get_variables(eom.natural_equation)
        # the Hopf variable (limit cycle frequency) does not contribute a term
        harm_vars = filter(x -> isequal(nat_var, x.natural_variable) && x.type !== "Hopf", eom.variables)
        ansatz = join([_show_ansatz(var) for var in harm_vars], " + ")
        output *= "\n" * string(nat_var) * " = " * ansatz
    end
    output
end


"""
$(TYPEDEF)

Holds a set of algebraic equations describing the steady state of a system.

# Fields
$(TYPEDFIELDS)

#  Constructors
```julia
Problem(eom::HarmonicEquation; Jacobian=true) # find and store the symbolic Jacobian
Problem(eom::HarmonicEquation; Jacobian="implicit") # ignore the Jacobian for now, compute implicitly later
Problem(eom::HarmonicEquation; Jacobian=J) # use J as the Jacobian (a function that takes a Dict)
Problem(eom::HarmonicEquation; Jacobian=false) # ignore the Jacobian
```
"""
mutable struct Problem
    "The harmonic variables to be solved for."
    variables::Vector{Num}
    "All symbols which are not the harmonic variables."
    parameters::Vector{Num}
    "The input object for HomotopyContinuation.jl solver methods."
    system::HomotopyContinuation.System
    "The Jacobian matrix (possibly symbolic).
    If `false`, the Jacobian is ignored (may be calculated implicitly after solving)."
    jacobian
    "The HarmonicEquation object used to generate this `Problem`."
    eom::HarmonicEquation

    Problem(variables,parameters,system,jacobian) = new(variables,parameters,system,jacobian) #incomplete initialization for user-defined symbolic systems
    Problem(variables,parameters,system,jacobian,eom) = new(variables,parameters,system,jacobian,eom)
end


function show(io::IO, p::Problem)
    println(io, length(p.system.expressions), " algebraic equations for steady states")
    println(io, "Variables: ", join(string.(p.variables), ", "))
    println(io, "Parameters: ", join(string.(p.parameters), ", "))
    println(io, "Symbolic Jacobian: ", !(p.jacobian==false))
end


"""
$(TYPEDEF)

Stores the steady states of a HarmonicEquation.

# Fields
$(TYPEDFIELDS)

"""
mutable struct Result
    "The variable values of steady-state solutions."
    solutions::Array{Vector{SteadyState}}
    "Values of all parameters for all solutions."
    swept_parameters::ParameterRange
    "The parameters fixed throughout the solutions."
    fixed_parameters::ParameterList
    "The `Problem` used to generate this."
    problem::Problem
    "Maps strings such as \"stable\", \"physical\" etc to arrays of values, classifying the solutions (see method `classify_solutions!`)."
    classes::Dict{String, Array}
    "The Jacobian with `fixed_parameters` already substituted. Accepts a dictionary specifying the solution.
    If problem.jacobian is a symbolic matrix, this holds a compiled function.
    If problem.jacobian was `false`, this holds a function that rearranges the equations to find J
    only after numerical values are inserted (preferable in cases where the symbolic J would be very large)."
    jacobian::Union{Function, Int64}

    Result(sol,swept, fixed, problem, classes, J) = new(sol, swept, fixed, problem, classes, J)
    Result(sol,swept, fixed, problem, classes) = new(sol, swept, fixed, problem, classes)
    Result(sol,swept, fixed, problem) = new(sol, swept, fixed, problem, Dict([]))
end


function show(io::IO, r::Result)
    println(io, "A steady state result for ", length(r.solutions), " parameter points")
    println(io, "\nSolution branches:   ", length(r.solutions[1]))
    println(io, "   of which real:    ", sum(any.(classify_branch(r, "physical"))))
    println(io, "   of which stable:  ", sum(any.(classify_branch(r, "stable"))))
    println(io, "\nClasses: ", join(keys(r.classes), ", "))
end


# overload to use [] for indexing
Base.getindex(r::Result, idx::Int...) = get_single_solution(r, idx)
Base.size(r::Result) = size(r.solutions)
