import HarmonicBalance.is_harmonic

@variables a,b,c,t,x(t),f, y(t)

@test is_harmonic(cos(f*t), t)
@test is_harmonic(1, t)
@test !is_harmonic(cos(f*t^2 + a), t)
@test !is_harmonic(a + t, t)

dEOM = DifferentialEquation([a + x, t^2 + cos(t)], [x, y])
@test !is_harmonic(dEOM, t)
