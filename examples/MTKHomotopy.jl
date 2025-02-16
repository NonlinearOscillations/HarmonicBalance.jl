using HarmonicBalance, ModelingToolkit, NonlinearSolve
using HarmonicBalance: OrderedDict
using NonlinearSolveHomotopyContinuation

begin
    rhs = function (u, p)
        return u * u - p[1] * u + p[2]
    end
    jac = function (u, p)
        return 2u - p[1]
    end
    fn = NonlinearFunction(rhs; jac)

    prob = NonlinearProblem(fn, 1.0, [5.0, 6.0])
    _alg = HomotopyContinuationJL{true}(; threading=false)
    sol = solve(prob, _alg)
end

@variables t x(t)
@parameters ω₀ γ λ F η α ω

natural_equation =
    d(d(x, t), t) +
    γ * d(x, t) +
    (ω₀^2 - λ * cos(2 * ω * t)) * x +
    α * x^3 +
    η * d(x, t) * x^2
forces = F * cos(ω * t)
diff_eq = DifferentialEquation(natural_equation + forces, x)

add_harmonic!(diff_eq, x, ω);
harmonic_eq = get_harmonic_equations(diff_eq)
eqs = getfield.(rearrange_standard(harmonic_eq).equations, :lhs) .~ 0

begin
    fixed = OrderedDict(
        ω₀ => 1.0, γ => 1e-2, λ => 5e-2, F => 1e-3, α => 1.0, η => 0.3, ω => 1.0
    )
    @mtkbuild ns = NonlinearSystem(eqs, get_variables(harmonic_eq), keys(fixed))

    fixed = OrderedDict(
        ω₀ => 1.0, γ => 1e-2, λ => 5e-2, F => 1e-3, α => 1.0, η => 0.3, ω => 1.0
    )
    fn = NonlinearFunction(ns; jac=true)
    hcp = NonlinearProblem(fn, [0.0, 0.0], fixed)
    solve(hcp)
    alg = HomotopyContinuationJL{true}(;)
    solve(hcp, alg)
end
begin
    @mtkbuild ns = NonlinearSystem(eqs, get_variables(harmonic_eq), harmonic_eq.parameters)

    fixed = (ω₀ => 1.0, γ => 1e-2, λ => 5e-2, F => 1e-3, α => 1.0, η => 0.3, ω => 1.0)
    hcp = HomotopyContinuationProblem(ns, [0.0, 0.0], fixed)
    solve(hcp)
    alg = HomotopyContinuationJL{true}(;)
    solve(hcp, alg)
end
