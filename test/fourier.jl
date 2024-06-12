import HarmonicBalance: fourier_cos_term, fourier_sin_term
import HarmonicBalance.Symbolics.expand

@variables f t θ a b

@test isequal(fourier_cos_term(cos(f * t)^2, f, t), 0)
@test isequal(fourier_sin_term(sin(f * t)^2, f, t), 0)

@test isequal(fourier_cos_term(cos(f * t)^2, 2 * f, t), 1//2)
@test isequal(fourier_sin_term(cos(f * t)^2, 2 * f, t), 0)
@test isequal(fourier_cos_term(sin(f * t)^2, 2 * f, t), -1//2)
@test isequal(fourier_sin_term(sin(f * t)^2, 2 * f, t), 0)

@test isequal(fourier_cos_term(cos(f * t), f, t), 1)
@test isequal(fourier_sin_term(sin(f * t), f, t), 1)

@test isequal(fourier_cos_term(cos(f * t + θ), f, t), cos(θ))
@test isequal(fourier_sin_term(cos(f * t + θ), f, t), -sin(θ))

to_check =
    (a * sin(f * t) + b * cos(f * t)) *
    (a * sin(2 * f * t) + b * cos(2 * f * t)) *
    (a * sin(3 * f * t) + b * cos(3 * f * t))
fourier_cos_term(to_check, 2 * f, t)
@test isequal(fourier_cos_term(to_check, 2 * f, t), expand(1//4 * (a^2 * b + b^3)))
@test isequal(fourier_cos_term(to_check, 4 * f, t), expand(1//4 * (a^2 * b + b^3)))
@test isequal(fourier_cos_term(to_check, 6 * f, t), expand(1//4 * (-3 * a^2 * b + b^3)))
@test isequal(fourier_sin_term(to_check, 2 * f, t), expand(1//4 * (a^3 + a * b^2)))
@test isequal(fourier_sin_term(to_check, 4 * f, t), expand(1//4 * (a^3 + a * b^2)))
@test isequal(fourier_sin_term(to_check, 6 * f, t), expand(1//4 * (-a^3 + 3 * a * b^2)))

# try something harder!
to_check = (a + b * cos(f * t + θ)^2)^3 * sin(f * t)
@test isequal(
    fourier_sin_term(to_check, f, t),
    expand(
        a^3 + a^2 * b * 3//2 + 9//8 * a * b^2 + 5//16 * b^3 -
        3//64 * b * (16 * a^2 + 16 * a * b + 5 * b^2) * cos(2 * θ),
    ),
)
@test isequal(
    fourier_cos_term(to_check, f, t),
    expand(-3//64 * b * (16 * a^2 + 16 * a * b + 5 * b^2) * sin(2 * θ)),
)

# FTing at zero : picking out constant terms
@test isequal(fourier_cos_term(cos(f * t), 0, t), 0)
@test isequal(fourier_cos_term(cos(f * t)^3 + 1, 0, t), 1)
@test isequal(fourier_cos_term(cos(f * t)^2 + 1, 0, t), 3//2)
@test isequal(fourier_cos_term((cos(f * t)^2 + cos(f * t))^3, 0, t), 23//16)
