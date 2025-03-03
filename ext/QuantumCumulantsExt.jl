module QuantumCumulantsExt

using QuantumCumulants: QuantumCumulants, average, MeanfieldEquations
using Symbolics: Symbolics, parse_expr_to_symbolic, Num, expand
using HarmonicBalance: HarmonicBalance
using OrderedCollections: OrderedDict

export Problem

function find_conjugate_pairs(ops)
    pairs = Tuple{Int,Int}[]
    for i in 1:length(ops)
        for j in i:length(ops)
            if isequal(ops[i], ops[j]')
                push!(pairs, (i, j))
            end
        end
    end
    return pairs
end
function get_filtered_equations(eqs::MeanfieldEquations)
    equations = getfield.(eqs.equations, :rhs)

    ops = getfield.(eqs.operator_equations, :lhs)
    conjugate_pairs = find_conjugate_pairs(ops)
    pairs = filter(pair -> !isequal(pair...), conjugate_pairs)

    foreach(pairs) do pair
        deleteat!(equations, pair[2])
    end
    return equations
end

function replace_to_reals(eqs::Vector{String}, ops_av::Vector{String}, conjugate_pairs)
    f(i) = iseven(i) ? "+" : "-"
    replacements = Pair{String,String}[]
    variables = String[]
    # prefactor = "√(2/(ħ))" # should contain frequency when we have multiple kinds of bosons

    foreach(conjugate_pairs) do pair
        is_real = isequal(pair...)
        op_pair = ops_av[pair[1]]
        op_name = replace(op_pair, "⟨" => "", "⟩" => "", "*" => "", "′" => "⁺")
        op_r = op_name * "ᵣ"
        op_i = op_name * "ᵢ"

        push!(variables, op_r)
        is_real || push!(variables, op_i)
        conj_replace = map(eachindex(pair)) do i
            op = ops_av[pair[i]] # ∨ f(i) contains "+" or "-"
            op => if is_real
                "((" * op_r * ") / √(2))"       #=prefactor * " " *=#
            else
                "((" * op_r * " $(f(i)) im*" * op_i * ") / √(2))"                #=prefactor *=#
            end
        end
        if isequal(pair...)
            push!(replacements, conj_replace[1])
        else
            push!(replacements, conj_replace...)
        end
    end

    return variables, replace.(eqs, replacements...)
end

function convert_to_symbolics(equations, ops_av, conjugate_pairs)
    eqs_str = string.(equations)
    ops_str = string.(ops_av)

    vars, eqs_re = replace_to_reals(eqs_str, ops_str, conjugate_pairs)
    eqs_expr = Meta.parse.(eqs_re)
    vars_expr = Meta.parse.(vars)
    eqs_symbolic = [parse_expr_to_symbolic(eq, @__MODULE__) for eq in eqs_expr]
    vars_symbolic = [parse_expr_to_symbolic(v, @__MODULE__) for v in vars_expr]
    return vars_symbolic, eqs_symbolic
end

function compute_real_equations(eqs::MeanfieldEquations)
    ops = union(eqs.operators, [c' for c in eqs.operators])
    ops_av = average.(ops)
    conjugate_pairs = find_conjugate_pairs(ops)

    equations = get_filtered_equations(eqs)
    vars, eqs_complex = convert_to_symbolics(equations, ops_av, conjugate_pairs)

    eqs_complex = expand.(eqs_complex)
    eqs_real = Num[]
    foreach(eqs_complex) do eq
        eq_re = √(2) * real(eq) # the factor of √(2) is needed to have Q to C limit.
        eq_im = -√(2) * imag(eq)
        iszero(eq_re) || push!(eqs_real, eq_re)
        iszero(eq_im) || push!(eqs_real, eq_im)
    end
    return vars, eqs_real
end # TEST: test if order of vars and eqs is correct

function HarmonicBalance.Problem(MFeqs::MeanfieldEquations, parameters, swept, fixed)
    vars, equations = compute_real_equations(MFeqs)

    return HarmonicBalance.Problem(
        Num.(equations), Num.(vars), Num.(parameters), OrderedDict(swept), OrderedDict(fixed)
    )
end

end
