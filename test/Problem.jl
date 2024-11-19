using HarmonicBalance
using Test

@testset "Problem without EOM" begin
    using HomotopyContinuation: System, Expression
    using Symbolics: Num, @variables
    F = System([Expression("x^2")])

    @variables x y
    prob = HarmonicBalance.Problem([x, y], Num[], F, [x y; x y])
    @test_throws UndefRefError prob.eom
end
