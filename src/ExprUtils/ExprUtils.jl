module ExprUtils

using DocStringExtensions
using OrderedCollections: OrderedDict

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
        simplify_fractions,
        substitute,
        term,
        expand,
        operation

    include("Symbolics_utils.jl")
    include("exponentials.jl")
    include("fourier.jl")
    include("drop_powers.jl")

end
