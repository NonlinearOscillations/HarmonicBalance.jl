using DispatchDoctor, JET
using HarmonicBalance
using Test

@testset "get_independent_variables" begin
    @variables ω1, ω2, ωₘ, t, ω, F, γ, λ, x(t), y(t)
    eqs = [d(x, t, 2) + (ω1^2 - λ * cos(ωₘ * t)) * x + γ * d(x, t)]

    diff_eq = DifferentialEquation(eqs, [x])
    @stable get_independent_variables(diff_eq)
end
@report_opt get_independent_variables(diff_eq)

@testset "get_all_terms" begin
    using HarmonicBalance.ExprUtils: get_all_terms, _get_all_terms
    using Symbolics
    @variables a, b, c
    @which get_all_terms()
    @code_warntype get_all_terms(a + b + c)
    @stable get_all_terms(a + b + c)
    expr = Symbolics.unwrap(a + b + c)
    @code_warntype _get_all_terms(expr)
end

@testset "get_independent" begin
    using HarmonicBalance.ExprUtils: get_independent
    @variables a, b, c, t

    @eqtest get_independent(a + b + c, t) == a + b + c
    @eqtest get_independent(a * b * c, t) == a * b * c
    @eqtest get_independent(a / b, t) == a / b
    @eqtest get_independent(a^2 + b^2 + c^2, t) == a^2 + b^2 + c^2
    @eqtest get_independent(a^2 / b^2, t) == a^2 / b^2
    @eqtest get_independent(2 * b^2, t) == 2 * b^2
    @eqtest get_independent(cos(t), t) == 0
    @eqtest get_independent(cos(t)^2 + 5, t) == 5
end
