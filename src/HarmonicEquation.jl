export get_independent_variables
export get_harmonic_equations
export is_rearranged
export slow_flow, slow_flow!
export _equations_without_brackets

show(eom::HarmonicEquation) = show_fields(eom)


"""
    harmonic_ansatz(eom::DifferentialEquation, time::Num; coordinates="Cartesian")

Expand each variable of `diff_eom` using the harmonics assigned to it with `time` as the time variable.
For each harmonic of each variable, an instance of `HarmonicVariable` (describing a pair of variables (u,v)) is automatically created and named. 

`coordinates` allows for using different coordinate systems (e.g. 'polars') - CURRENTLY INACTIVE
"""
harmonic_ansatz(eom::DifferentialEquation, time::Num; coordinates="Cartesian") = harmonic_ansatz(eom, time, Dict([[var, coordinates] for var in get_variables(eom)]))


function harmonic_ansatz(diff_eom::DifferentialEquation, time::Num, trafo_types::Dict{Num,String})
    !is_harmonic(diff_eom, time) && error("The differential equation is not harmonic in ", time, " !")
    old_eqs = [diff_eom.equations[var] for var in get_variables(diff_eom)]
    eqs, vars = deepcopy(old_eqs), []
    rules = Dict()

    idx_counter = 1 # keep count to label new variables
    for var in get_variables(diff_eom)
        to_substitute = 0 # holds all the subtitution rules
        for ω in diff_eom.harmonics[var] # 2 new variables are created for each ω of each variable
            these_names = coordinate_names[trafo_types[var]] # stores two strings denoting the variables as either "a", "ϕ", "u", "v"
            unique_symbols = [  these_names[1] * string(idx_counter) , these_names[2] * string(idx_counter) ]
            repl_rule, new_vars = _rotate_variable(var, ω, time, trafo_types[var], new_symbols=unique_symbols) # create the new variables
            to_substitute += repl_rule
            append!(vars, [new_vars]) 
            idx_counter += 1
        end
        rules[var] = to_substitute
    end
    eqs = substitute_all(eqs, rules)
    HarmonicEquation(eqs, Vector{HarmonicVariable}(vars), diff_eom)
end 


function slow_flow!(eom::HarmonicEquation; fast_time::Num, slow_time::Num, degree=2)
    eom.equations = expand_derivatives.(eom.equations) # expand all the derivatives

    # fast_time => slow_time for derivatives up to degree-1
    replace0 = [var => substitute_all(var, fast_time => slow_time) for var in get_variables(eom)]
    replace_degrees = [Differential(fast_time)^deg => Differential(slow_time)^deg for deg in 1:degree-1]
    replace = flatten([replace0, replace_degrees])

    # degree derivatives are removed
    drop = [d(var, fast_time, degree) => 0 for var in get_variables(eom)]

    eom.equations = substitute_all(substitute_all(eom.equations, drop), replace)
    eom.variables = substitute_all(eom.variables, replace)
end


"""
    slow_flow(eom::HarmonicEquation; fast_time::Num, slow_time::Num, degree=2)

Removes all derivatives w.r.t `fast_time` (and their products) in `eom` of power `degree`.
In the remaining derivatives, `fast_time` is replaced by `slow_time`.
"""
function slow_flow(eom::HarmonicEquation; fast_time::Num, slow_time::Num, degree=2)::HarmonicEquation
    new_eq = deepcopy(eom)
    slow_flow!(new_eq, fast_time=fast_time, slow_time=slow_time, degree=degree)
    new_eq
end


#   Drop powers of `var` of degree >= `deg` from the equation set in `eom`.
function drop_powers(eom::HarmonicEquation, var, deg::Integer)
    new_eom = deepcopy(eom)
    new_eom.equations = drop_powers(eom.equations, var, deg)
    new_eom
end


"Rearrange an equation system such that the field equations is equal to the vector specified in new_lhs"
function rearrange!(eom::HarmonicEquation, new_rhs::Vector{Num})
    soln = Symbolics.solve_for(eom.equations, new_rhs, simplify=false,check=true)
    eom.equations = soln .~ new_rhs
    return nothing
end


function rearrange(eom::HarmonicEquation, new_rhs::Vector{Num})
    new_eom = deepcopy(eom)
    rearrange!(new_eom, new_rhs)
    return new_eom
end


"""
$(TYPEDSIGNATURES)
Check if `eom` is rearranged to the standard form, such that the derivatives of the variables are on one side.
"""
function is_rearranged(eom::HarmonicEquation)
    tvar = get_independent_variables(eom)[1]
    isequal(getfield.(eom.equations, :rhs), d(get_variables(eom), tvar))
end


"""
$(TYPEDSIGNATURES)
Rearrange `eom` to the standard form, such that the derivatives of the variables are on one side.
"""
function rearrange_standard(eom::HarmonicEquation)
    tvar = get_independent_variables(eom)[1]
    dvars =  d(get_variables(eom), tvar)
    rearrange(eom, dvars)
end


"""
$(TYPEDSIGNATURES)
Get the internal symbols of the independent variables of `eom`. 
"""
function get_variables(eom::HarmonicEquation)
    return flatten(get_variables.(eom.variables))
end


"Get the parameters (not time nor variables) of a HarmonicEquation"
function _parameters(eom::HarmonicEquation)
    all_symbols = flatten([cat(get_variables(eq.lhs), get_variables(eq.rhs), dims=1) for eq in eom.equations])
    # subtract the set of independent variables (i.e., time) from all free symbols
    setdiff(all_symbols, get_variables(eom), get_independent_variables(eom)) 
end


###
# Extending Symbolics.jl's simplify and substitute
###

"Apply `rules` to both `equations` and `variables` field of `eom`"
function substitute_all!(eom::HarmonicEquation, rules::Dict)
    eom.equations = replace(eom.equations, rules)
    eom.variables = substitute_all(eom.variables, rules)
end


"Simplify the equations in HarmonicEquation."
function simplify!(eom::HarmonicEquation)
    eom.equations = [simplify(eq) for eq in eom.equations]
end



"""
$(TYPEDSIGNATURES)
Return the independent variables (typically time) of `eom`.
"""
function get_independent_variables(eom::HarmonicEquation)::Vector{Num}
    dynamic_vars = flatten(getfield.(eom.variables, Symbol("symbols")))
    return flatten(unique([SymbolicUtils.arguments(var.val) for var in dynamic_vars]))
end


"""
$(TYPEDSIGNATURES)
Extract the Fourier components of `eom` corresponding to the harmonics specified in `eom.variables`.
For each harmonic of each variable, 2 equations are generated (cos and sin Fourier coefficients).
`time` does not appear in the resulting equations anymore.

Underlying assumption: all time-dependences are harmonic.
"""
function fourier_transform(eom::HarmonicEquation, time::Num)
    new_eom = deepcopy(eom)
    fourier_transform!(new_eom, time)
    return new_eom
end


function fourier_transform!(eom::HarmonicEquation, time::Num)
    avg_eqs = Vector{Equation}(undef, 2*length(eom.variables))

    natural_variables = flatten([[x,x] for x in getfield.(eom.variables, :natural_variable)])
    equation_assignments = flatten([findall(x -> isequal(x, var), collect(keys(eom.natural_equation.equations))) for var in natural_variables])
    
    for i in 1:length(avg_eqs)
        j = Int(ceil(i/2))
        ω = eom.variables[j].ω
        eq = eom.equations[equation_assignments[i]]

        if isodd(i)
            avg_eqs[i] = fourier_cos_term(eq, ω, time)
        else
            avg_eqs[i] = fourier_sin_term(eq, ω, time)
        end
    end
    
    eom.equations = avg_eqs
end


"""
    get_harmonic_equations(diff_eom::DifferentialEquation; fast_time=nothing, slow_time=nothing)

Apply the harmonic ansatz, followed by the slow-flow, Fourier transform and dropping 
higher-order derivatives to obtain
a set of ODEs (the harmonic equations) governing the harmonics of `diff_eom`.

The harmonics evolve in `slow_time`, the oscillating terms themselves in `fast_time`.
If no input is used, a variable T is defined for `slow_time` and `fast_time` is taken as the independent variable
of `diff_eom`.

By default, all products of order > 1 of `slow_time`-derivatives are dropped,
which means the equations are linear in the time-derivatives.

# Example
```julia-repl
julia> @variables t, x(t), ω0, ω, F;

# enter the simple harmonic oscillator
julia> diff_eom = DifferentialEquation( d(x,t,2) + ω0^2 * x ~ F *cos(ω*t), x);

# expand x in the harmonic ω
julia> add_harmonic!(diff_eom, x, ω);

# get equations for the harmonics evolving in the slow time T
julia> harmonic_eom = get_harmonic_equations(diff_eom)

A set of 2 harmonic equations
Variables: u1(T), v1(T)
Parameters: ω0, ω, F

Harmonic ansatz: 
x(t) = u1*cos(ωt) + v1*sin(ωt)

Harmonic equations:

(ω0^2)*u1(T) + (2//1)*ω*Differential(T)(v1(T)) - (ω^2)*u1(T) ~ F

(ω0^2)*v1(T) - (ω^2)*v1(T) - (2//1)*ω*Differential(T)(u1(T)) ~ 0
```

"""
function get_harmonic_equations(diff_eom::DifferentialEquation; fast_time=nothing, slow_time=nothing)

    slow_time = isnothing(slow_time) ? (@variables T; T) : slow_time
    fast_time = isnothing(fast_time) ? get_independent_variables(diff_eom)[1] : fast_time

    all(isempty.(values(diff_eom.harmonics))) && error("No harmonics specified!")
    eom = harmonic_ansatz(diff_eom, fast_time); # substitute trig functions into the differential equation
    eom = slow_flow(eom, fast_time=fast_time, slow_time=slow_time); # drop 2nd order time derivatives
    fourier_transform!(eom, fast_time); # perform averaging over the frequencies originally specified in dEOM
    ft_eom_simplified = drop_powers(eom, d(get_variables(eom), slow_time), 2); # drop higher powers of the first-order derivatives
    return ft_eom_simplified 
end


"Rearrange `eq` to have zero on the right-hand-side."
_set_zero_rhs(eq::Equation) = eq.lhs - eq.rhs ~ 0
_set_zero_rhs(eqs::Vector{Equation}) = [_set_zero_rhs(eq) for eq in eqs]


"Returns the equation system in `eom`, dropping all argument brackets (i.e., u(T) becomes u)."
function _equations_without_brackets(eom::HarmonicEquation)
    variable_rules = [var => declare_variable(var_name(var)) for var in get_variables(eom)]
    equations_lhs = Num.(getfield.(eom.equations, :lhs)  - getfield.(eom.equations, :rhs))
    substitute_all(equations_lhs, variable_rules)
end


function drop_powers(eom::HarmonicEquation, terms, deg::Int)
    new_eom = deepcopy(eom)
    new_eom.equations = drop_powers(eom.equations, terms, deg)
    new_eom
end