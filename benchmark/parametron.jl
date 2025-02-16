using HarmonicBalance, Plots
using BenchmarkTools
using HomotopyContinuation

@variables ω₀ γ λ α ω t x(t)

natural_equation = d(d(x, t), t) + γ * d(x, t) + (ω₀^2 - λ * cos(2 * ω * t)) * x + α * x^3
diff_eq = DifferentialEquation(natural_equation, x)

add_harmonic!(diff_eq, x, ω);

harmonic_eq = get_harmonic_equations(diff_eq)
# harmonic_eq = HarmonicBalance.harmonic_eqlem(harmonic_eq)

fixed = (ω₀ => 1.0, γ => 0.002, α => 1.0)
varied = (ω => range(0.99, 1.01, 100), λ => range(1e-6, 0.03, 100))

track_options = HomotopyContinuation.TrackerOptions(;
    parameters=:conservative,
    automatic_differentiation=3,
    max_steps=20_000,
    min_step_size=1e-48,
)
end_options = HomotopyContinuation.EndgameOptions(; refine_steps=10)
method = WarmUp(; compile=true, tracker_options=track_options, endgame_options=end_options)
result_2D = get_steady_states(harmonic_eq, method, varied, fixed; show_progress=false)
plot_phase_diagram(result_2D; class="stable")

# loose some solutions
track_options = HomotopyContinuation.TrackerOptions(;
    parameters=:fast, extended_precision=false
)
end_options = HomotopyContinuation.EndgameOptions(;
    endgame_start=0.0, only_nonsingular=true
)
method = WarmUp(; compile=true, tracker_options=track_options, endgame_options=end_options)
result_2D = get_steady_states(harmonic_eq, method, varied, fixed; show_progress=false)
plot_phase_diagram(result_2D; class="stable")

@btime $get_steady_states($harmonic_eq, $WarmUp(), $varied, $fixed; show_progress=false)
# 1.029 s (4495544 allocations: 445.47 MiB)

@btime $get_steady_states(
    $harmonic_eq, $WarmUp(; compile=true), $varied, $fixed; show_progress=false
) # 953.323 ms (4490408 allocations: 442.09 MiB)

track_options = HomotopyContinuation.TrackerOptions(; parameters=:fast)
method = WarmUp(; compile=true, tracker_options=track_options)
@btime $get_steady_states($harmonic_eq, $method, $varied, $fixed; show_progress=false) #  747.007 ms (4460916 allocations: 438.87 MiB)

track_options = HomotopyContinuation.TrackerOptions(; parameters=:fast)
end_options = HomotopyContinuation.EndgameOptions(; endgame_start=0.0)
method = WarmUp(; compile=true, tracker_options=track_options, endgame_options=end_options)
@btime $get_steady_states($harmonic_eq, $method, $varied, $fixed; show_progress=false) # 734.310 ms (4460118 allocations: 438.81 MiB)

track_options = HomotopyContinuation.TrackerOptions(;
    parameters=:fast, extended_precision=false
)
end_options = HomotopyContinuation.EndgameOptions(; endgame_start=0.0)
method = WarmUp(; compile=true, tracker_options=track_options, endgame_options=end_options)
@btime $get_steady_states($harmonic_eq, $method, $varied, $fixed; show_progress=false) # 734.310 ms (4460118 allocations: 438.81 MiB)
