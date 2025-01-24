module ExprUtils

using DocStringExtensions
using OrderedCollections: OrderedDict
using LinearAlgebra: LinearAlgebra

using SymbolicUtils:
    SymbolicUtils,
    Postwalk,
    Sym,
    BasicSymbolic,
    isterm,
    ispow,
    isadd,
    isdiv,
    ismul,
    add_with_div,
    frac_maketerm,
    @compactified,
    issym

using Symbolics:
    Symbolics,
    Num,
    unwrap,
    wrap,
    get_variables,
    Equation,
    Differential,
    @variables,
    arguments,
    substitute,
    term,
    expand,
    operation

include("Symbolics_utils.jl")
include("exponentials.jl")
include("fourier.jl")
include("drop_powers.jl")

macro eqtest(expr)
    @assert expr.head == :call && expr.args[1] in [:(==), :(!=)]
    return esc(
        if expr.args[1] == :(==)
            :(@test isequal($(expr.args[2]), $(expr.args[3])))
        else
            :(@test !isequal($(expr.args[2]), $(expr.args[3])))
        end,
    )
end

macro eqsym(expr)
    @assert expr.head == :call && expr.args[1] in [:(==), :(!=)]
    return esc(
        if expr.args[1] == :(==)
            :(isequal($(expr.args[2]), $(expr.args[3])))
        else
            :(!isequal($(expr.args[2]), $(expr.args[3])))
        end,
    )
end

is_identity(A::Matrix{Num}) = (@eqsym A == Matrix{Num}(LinearAlgebra.I, size(A)...))

end
