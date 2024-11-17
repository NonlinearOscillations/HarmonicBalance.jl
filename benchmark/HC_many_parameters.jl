using HomotopyContinuation
@var u1, v1, ω, α, γ, λ, ω0

eqs = [
    -u1 * ω^2 +
    u1 * ω0^2 +
    (3 / 4) * u1^3 * α +
    (3 / 4) * u1 * v1^2 * α +
    (-1 / 2) * u1 * λ * ω0^2 +
    v1 * γ * ω,
    -v1 * ω^2 + v1 * ω0^2 + (3 / 4) * v1^3 * α - u1 * γ * ω +
    (3 / 4) * u1^2 * v1 * α +
    (1 / 2) * v1 * λ * ω0^2,
]

F = System(eqs; parameters=[ω, α, γ, λ, ω0], variables=[u1, v1])

input_array = [
    [0.9, 1.0, 0.00, 0.05, 1.0],
    [0.93, 1.0, 0.00, 0.05, 1.0],
    [0.97, 1.0, 0.0, 0.05, 1.0],
    [1.0, 1.0, 0.00, 0.05, 1.0],
    [1.03, 1.0, 0.00, 0.05, 1.0],
    [1.06, 1.0, 0.00, 0.05, 1.0],
]
generic_parameters = randn(ComplexF64, 5)

R0 = solve(F; target_parameters=generic_parameters, threading=true)
R1 = solve(
    F,
    solutions(R0);
    start_parameters=generic_parameters,
    target_parameters=input_array,
    threading=true,
    only_non_zero=true,
)
solutions(R1[4][1])
solutions(R0)
R1[4]
