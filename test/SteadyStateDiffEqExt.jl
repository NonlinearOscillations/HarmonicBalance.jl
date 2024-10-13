
using HarmonicBalance, ModelingToolkit, SteadyStateDiffEq, OrdinaryDiffEqRosenbrock, LinearAlgebra, NonlinearSolve

@testset "Steady state sweeps" begin
    @testset "one variable ODE" begin
        @independent_variables t
        @variables v(t) = 0
        @parameters g = 9.8 k = 0.2
        D = Differential(t)
        eqs = [D(v) ~ g - k * v]
        @named model = ODESystem(eqs, t)

        model = structural_simplify(model)

        prob_ss = SteadyStateProblem{false}(model, [], []; jac=true)
        prob_np = NonlinearProblem(prob_ss)

        @test prob_np.p isa MTKParameters
        @test prob_np.p.tunable isa Vector{Float64}

        varied = 2 => range(0, 1, 100)

        swept = steady_state_sweep(
            prob_np, prob_ss, NewtonRaphson(), DynamicSS(Rodas5()); varied=varied
        )
    end

    @testset "two variable ODE (duffing)" begin
        @independent_variables t
        @variables u1(t) v1(t)
        @parameters α ω ω0 F γ

        eqs = [
            Differential(t)(u1) ~
                (
                    F * γ - u1 * γ * (ω^2) - u1 * γ * (ω0^2) - v1 * (γ^2) * ω -
                    2 * v1 * (ω^3) + 2 * v1 * ω * (ω0^2) - (3//4) * (u1^3) * α * γ +
                    (3//2) * (u1^2) * v1 * α * ω - (3//4) * u1 * (v1^2) * α * γ +
                    (3//2) * (v1^3) * α * ω
                ) / (γ^2 + (4) * (ω^2)),
            Differential(t)(v1) ~
                (
                    -2 * F * ω - u1 * (γ^2) * ω - (2) * u1 * (ω^3) +
                    2 * u1 * ω * (ω0^2) +
                    v1 * γ * (ω^2) +
                    v1 * γ * (ω0^2) +
                    (3//2) * (u1^3) * α * ω +
                    (3//4) * (u1^2) * v1 * α * γ +
                    (3//2) * u1 * (v1^2) * α * ω +
                    (3//4) * (v1^3) * α * γ
                ) / (-(γ^2) - (4) * (ω^2)),
        ]

        @named model = ODESystem(eqs, t, [u1, v1], [α, ω, ω0, F, γ])
        model = structural_simplify(model)

        param =  [α, ω, ω0, F, γ] .=> [1.0, 1.2, 1.0, 0.01, 0.01]
        x0 = [1.0, 0.0]
        prob_ss = SteadyStateProblem{true}(model, x0, param; jac=true)
        prob_np = NonlinearProblem(prob_ss)

        ω_span = (0.9, 1.5)
        ω_range = range(ω_span..., 100)
        varied_idx = findfirst(x -> x == 1.2, prob_ss.p.tunable)
        varied = varied_idx => ω_range
        swept = steady_state_sweep(
            prob_np, prob_ss, NewtonRaphson(), DynamicSS(Rodas5()); varied=varied
        )

        @test length(swept) == 100
        @test length(swept[1]) == 2
        @test norm.(swept) isa Array{Float64,1}
        # using Plots, LinearAlgebra; plot(norm.(swept))

        function has_discontinuity(v::Vector{Float64})
            threshold = 1.0e-1  # Define a threshold for the discontinuity
            for i in 2:length(v)
                abs(v[i] - v[i - 1]) > threshold && return true
            end
            return false
        end

        @test has_discontinuity(norm.(swept))
    end
end # Steady state sweeps
