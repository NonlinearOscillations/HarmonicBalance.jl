@testset "Symbolics arrays" begin
    using Symbolics, HarmonicBalance
    @variables x[1:2]
    vars = Symbolics.scalarize(x)
    HarmonicBalance.var_name(vars[1])
end
