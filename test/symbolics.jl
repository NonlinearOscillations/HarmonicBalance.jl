using Test
using Symbolics
using HarmonicBalance
using SymbolicUtils: Fixpoint, Prewalk, PassThrough

macro eqtest(expr)
    @assert expr.head == :call && expr.args[1] in [:(==), :(!=)]
    if expr.args[1] == :(==)
        :(@test isequal($(expr.args[2]), $(expr.args[3])))
    else
        :(@test !isequal($(expr.args[2]), $(expr.args[3])))
    end |> esc
end

@testset "exp(x)^n => exp(x*n)" begin
    using HarmonicBalance: expand_all
    @variables a n

    @eqtest simplify(exp(a)^3) == exp(3 * a)
    @eqtest simplify(exp(a)^n) == exp(n * a)
    @eqtest expand_all(exp(a)^3) == exp(3 * a)
    @eqtest expand_all(exp(a)^3) == exp(3 * a)
    @eqtest expand_all(im * exp(a)^5) == im * exp(5 * a)
end

@testset "exp(a)*exp(b) => exp(a+b)" begin
    using HarmonicBalance: simplify_exp_products
    @variables a b

    @eqtest simplify_exp_products(exp(a) * exp(b)) == exp(a + b)
    @eqtest simplify_exp_products(exp(3a) * exp(4b)) == exp(3a + 4b)
    @eqtest simplify_exp_products(im * exp(3a) * exp(4b)) == im * exp(3a+4b)
end


@testset "simplify_complex" begin
    using HarmonicBalance: simplify_complex, is_not_complex
    @variables a b

    z = Complex{Num}(a)
    @test !is_not_complex(a)
    @test is_not_complex(z)
    @eqtest simplify_complex(z) == a
    @test a isa Num

    z = Complex{Num}(1+0*im)
    @test simplify_complex(z) isa Num
end

@testset "euler" begin
    @variables a b
    @eqtest cos(a) + im*sin(a) == exp(im*a)
    @eqtest exp(a)*cos(b) + im*sin(b)*exp(a) == exp(a + im*b)
end

@testset "powers" begin
    using HarmonicBalance: drop_powers, max_power
    using HarmonicBalance.Symbolics: expand

    @variables a, b, c

    @eqtest max_power(a^2 + b, a) == 2
    @eqtest max_power(a * ((a + b)^4)^2 + a, a) == 9

    @eqtest drop_powers(a^2 + b, a, 1) == b
    @eqtest drop_powers((a + b)^2, a, 1) == b^2
    @eqtest drop_powers((a + b)^2, [a, b], 1) == 0
    @eqtest drop_powers((a + b)^3 + (a + b)^5, [a, b], 4) == expand((a + b)^3)
end

@testset "trig_to_exp" begin
    using HarmonicBalance: expand_all, trig_to_exp, unwrap
    using SymbolicUtils.Code: toexpr
    @variables f t
    cos_euler(x) = (exp(im*x)+exp(-im*x))/2
    sin_euler(x) = (exp(im*x)-exp(-im*x))/(2*im)
    @eqtest string(trig_to_exp(cos(f*t))) == "(1//2)*(exp((0 - 1im)*f*t) + exp((0 + 1im)*f*t))"
    @eqtest string(trig_to_exp(sin(f*t))) == "(0.0 - 0.5im)*(-exp((0 - 1im)*f*t) + exp((0 + 1im)*f*t))"
    @eqtest string(trig_to_exp(sin(f*t)^2 + cos(f*t))) == "(1//2)*(exp((0 - 1im)*f*t) + exp((0 + 1im)*f*t)) + (-0.25 - 0.0im)((-exp((0 - 1im)*f*t) + exp((0 + 1im)*f*t))^2)"
end

@testset "harmonic" begin
    using HarmonicBalance: is_harmonic

    @variables a, b, c, t, x(t), f, y(t)

    @test is_harmonic(cos(f * t), t)
    @test is_harmonic(1, t)
    @test !is_harmonic(cos(f * t^2 + a), t)
    @test !is_harmonic(a + t, t)

    dEOM = DifferentialEquation([a + x, t^2 + cos(t)], [x, y])
    @test !is_harmonic(dEOM, t)
end

# @testset "fourier" begin
#     using HarmonicBalance: fourier_cos_term, fourier_sin_term
#     using HarmonicBalance.Symbolics: expand
#     using HarmonicBalance: add_div, trig_to_exp

#     @variables f t θ a b
#     x = cos(f * t)^2
#     x = add_div(x) # a/b + c/d = (ad + bc)/bd
#     x = expand(x) # open all brackets
#     x = trig_to_exp(x)
#     x = simplify(x)

#     @test isequal(fourier_cos_term(cos(f * t)^2, f, t), 0)
#     @test isequal(fourier_sin_term(sin(f * t)^2, f, t), 0)

#     @test isequal(fourier_cos_term(cos(f * t)^2, 2 * f, t), 1//2)
#     @test isequal(fourier_sin_term(cos(f * t)^2, 2 * f, t), 0)
#     @test isequal(fourier_cos_term(sin(f * t)^2, 2 * f, t), -1//2)
#     @test isequal(fourier_sin_term(sin(f * t)^2, 2 * f, t), 0)

#     @test isequal(fourier_cos_term(cos(f * t), f, t), 1)
#     @test isequal(fourier_sin_term(sin(f * t), f, t), 1)

#     @test isequal(fourier_cos_term(cos(f * t + θ), f, t), cos(θ))
#     @test isequal(fourier_sin_term(cos(f * t + θ), f, t), -sin(θ))

#     term =
#         (a * sin(f * t) + b * cos(f * t)) *
#         (a * sin(2 * f * t) + b * cos(2 * f * t)) *
#         (a * sin(3 * f * t) + b * cos(3 * f * t))
#     fourier_cos_term(term, 2 * f, t)
#     @test isequal(fourier_cos_term(term, 2 * f, t), expand(1//4 * (a^2 * b + b^3)))
#     @test isequal(fourier_cos_term(term, 4 * f, t), expand(1//4 * (a^2 * b + b^3)))
#     @test isequal(fourier_cos_term(term, 6 * f, t), expand(1//4 * (-3 * a^2 * b + b^3)))
#     @test isequal(fourier_sin_term(term, 2 * f, t), expand(1//4 * (a^3 + a * b^2)))
#     @test isequal(fourier_sin_term(term, 4 * f, t), expand(1//4 * (a^3 + a * b^2)))
#     @test isequal(fourier_sin_term(term, 6 * f, t), expand(1//4 * (-a^3 + 3 * a * b^2)))

#     # try something harder!
#     term = (a + b * cos(f * t + θ)^2)^3 * sin(f * t)
#     @test isequal(
#         fourier_sin_term(term, f, t),
#         expand(
#             a^3 + a^2 * b * 3//2 + 9//8 * a * b^2 + 5//16 * b^3 -
#             3//64 * b * (16 * a^2 + 16 * a * b + 5 * b^2) * cos(2 * θ),
#         ),
#     )
#     @test isequal(
#         fourier_cos_term(term, f, t),
#         expand(-3//64 * b * (16 * a^2 + 16 * a * b + 5 * b^2) * sin(2 * θ)),
#     )

#     # FTing at zero : picking out constant terms
#     @test isequal(fourier_cos_term(cos(f * t), 0, t), 0)
#     @test isequal(fourier_cos_term(cos(f * t)^3 + 1, 0, t), 1)
#     @test isequal(fourier_cos_term(cos(f * t)^2 + 1, 0, t), 3//2)
#     @test isequal(fourier_cos_term((cos(f * t)^2 + cos(f * t))^3, 0, t), 23//16)
# # end
