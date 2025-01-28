using HarmonicBalance
using HarmonicBalance: OrderedDict
using Test

@testset "Problem without EOM" begin
    using HomotopyContinuation: System, Expression
    using Symbolics: Num, @variables
    F = System([Expression("x^2")])

    @variables x y
    prob = HarmonicBalance.Problem(
        [x, y],
        Num[],
        OrderedDict{Num,Vector{Float64}}(),
        OrderedDict{Num,Float64}(),
        F,
        HarmonicBalance.JacobianFunction(ComplexF64)(x -> x),
    )
    @test_throws UndefRefError prob.eom

    prob = HarmonicBalance.Problem(
        [x, y], Num[], OrderedDict{Num,Vector{Float64}}(), OrderedDict{Num,Float64}(), F
    )
    @test_throws UndefRefError prob.jacobian
end
