using HarmonicBalance;
HB = HarmonicBalance;

@variables t T x(t) y(t) # symbolic variables
@variables ω ω0 γ F α λ ψ θ η

eq = [
    d(d(x, t), t) +
    γ * d(x, t) +
    ω0^2 * (1 - λ * cos(2 * ω * t + ψ)) * x +
    α * x^3 +
    η * d(x, t) * x^2 ~ F * cos(ω * t + θ),
]

diff_eom = DifferentialEquation(eq, [x])

add_harmonic!(diff_eom, x, ω) # x will rotate at ω

harmonic_eq1 = get_krylov_equations(diff_eom; order=1)
harmonic_eq2 = get_krylov_equations(diff_eom; order=2)

@testset "show method" begin
    print_variable = HB._show_ansatz(harmonic_eq1)
    @test !(occursin("xˍt(t)", print_variable))
end

fixed = (ω0 => 1.0, γ => 0.005, α => 1.0, η => 0, F => 0.0, ψ => 0.0, θ => 0.0)
varied = (ω => range(0.99, 1.01, 5), λ => range(1e-6, 0.05, 5))

res1 = get_steady_states(harmonic_eq1, varied, fixed; show_progress=false, seed=SEED);
res2 = get_steady_states(harmonic_eq2, varied, fixed; show_progress=false, seed=SEED);

@testset "not damped" begin
    using HarmonicBalance.ExprUtils: expand_fraction

    @testset "single resonator" begin
        @variables t x(t) y(t) ω0 ω F α # symbolic variables
        eq1 = d(d(x, t), t) + ω0^2 * x + α * x^3 ~ F * cos(ω * t)
        EOM = DifferentialEquation(eq1, x)
        add_harmonic!(EOM, x, ω)
        krylov_eq = get_krylov_equations(EOM; order=1)
        harmonic_eq = get_harmonic_equations(EOM)
        rearranged = HarmonicBalance.rearrange_standard(harmonic_eq)

        @testset for i in 1:2
            eqk = expand_fraction(krylov_eq.equations[i].lhs)
            eqh = expand_fraction(rearranged.equations[i].lhs)
            @variables T u1(T) v1(T) ω0 ω F α # symbolic variables
            subs = Dict(u1 => 1, v1 => 1, α => 1, F => 1, ω0 => 1, ω => 1)
            solk = substitute(eqk, subs)
            solh = substitute(eqh, subs)
            @test Float64(solk + solh) == 0.0
            # ^ different ansatz
        end
    end

    @testset "multiple resonators" begin
        @variables t x(t) y(t) ω0 ω F α J # symbolic variables
        eq1 = d(d(x, t), t) + ω0^2 * x + α * x^3 ~ F * cos(ω * t) + J * y
        eq2 = d(d(y, t), t) + ω0^2 * y + α * y^3 ~ F * cos(ω * t) + J * x
        EOM = DifferentialEquation([eq1, eq2], [x, y])
        add_harmonic!(EOM, x, ω)
        add_harmonic!(EOM, y, ω)
        krylov_eq = get_krylov_equations(EOM; order=1)
        harmonic_eq = get_harmonic_equations(EOM)
        rearranged = HarmonicBalance.rearrange_standard(harmonic_eq)

        @testset for i in 1:4
            eqk = expand_fraction(krylov_eq.equations[i].lhs)
            eqh = expand_fraction(rearranged.equations[i].lhs)
            @variables T u1(T) v1(T) ω0 ω F α # symbolic variables
            subs = Dict([u1, v1, α, F, ω0, ω, J] .=> rand(7))
            solk = substitute(eqk, subs)
            solh = substitute(eqh, subs)
            @test Float64(solk + solh) == 0.0
            # ^ different ansatz
        end
    end

    @testset "two resonators at different frequencies" begin
        @variables ω₁, ω₂, t, ω, F, γ₁, γ₂
        @variables α₁, J₂, J₁, α₂, η₁, η₂
        @variables x(t), y(t)
        # eq1 = d(d(x, t), t) + ω₁^2 * x + α₁ * x^3 + 3 * J₁ * x^2 * y + J₂ * y^2 * x
        # eq2 = d(d(y, t), t) + ω₂^2 * y + α₂ * y^3 + J₁ * x^3 + J₂ * x^2 * y
        eq1 = d(d(x, t), t) + ω₁^2 * x + α₁ * x^3 + 3 * J₁ * x^2 * y + J₂ * y^2 * x
        eq2 = d(d(y, t), t) + ω₂^2 * y + α₂ * y^3 + J₁ * x^3 + J₂ * x^2 * y
        forces = [F * cos(ω * t), 0]
        dEOM_temp = DifferentialEquation([eq1, eq2] - forces, [x, y])

        add_harmonic!(dEOM_temp, x, ω) # x will rotate at ω
        add_harmonic!(dEOM_temp, y, 3 * ω) # y will rotate at 3*ω

        krylov_eq = get_krylov_equations(dEOM_temp; order=1)
        harmonic_eq = get_harmonic_equations(dEOM_temp)
        rearranged = HarmonicBalance.rearrange_standard(harmonic_eq)

        @testset for i in 1:4
            eqk = expand_fraction(krylov_eq.equations[1].lhs)
            eqh = expand_fraction(rearranged.equations[1].lhs)
            @variables T u1(T) v1(T) u2(T) v2(T)
            subs = Dict([u1, v1, u2, v2, ω₁, ω₂, ω, F, J₂, J₁, α₁, α₂] .=> rand(12))
            solk = substitute(eqk, subs)
            solh = substitute(eqh, subs)
            @test Float64(solk + solh) == 0.0 broken = true
            # ^ different ansatz
        end
    end
end
