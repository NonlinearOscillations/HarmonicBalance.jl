using HarmonicBalance

# define equation of motion
@variables ω1, ω2, t, ω, F, γ, α1, α2, k, x(t), y(t);
rhs = [
    d(x, t, 2) + ω1^2 * x + γ * d(x, t) + α1 * x^3 - k * y,
    d(d(y, t), t) + ω2^2 * y + γ * d(y, t) + α2 * y^3 - k * x,
]
eqs = rhs .~ [F * cos(ω * t), 0]

@test_throws ArgumentError DifferentialEquation(rhs ~ [F * cos(ω * t), 0], [x, y])
