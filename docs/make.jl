push!(LOAD_PATH, "../src/")

using Documenter
using HarmonicBalance

makedocs(
	sitename="HarmonicBalance.jl",
	format = Documenter.HTML(mathengine=MathJax()),
	pages = [
		"Background" => Any[
			"background/harmonic_balance.md"
			"background/stability_response.md"
		],
		"Examples" => Any[
			"examples/single_Duffing.md"
			"examples/linear_response.md"
			"examples/time_dependent.md"
			"examples/single_parametron_1D.md"
			"examples/single_parametron_2D.md"
			"examples/single_parametron_time_dep.md"
			],
		"Manual" => Any[
			"manual/entering_eom.md"
			"manual/extracting_harmonics.md"
			"manual/solving_harmonics.md"
			"manual/plotting.md"
			"manual/time_dependent.md"
			"manual/linear_response.md"
			"manual/saving.md"
			]
		]
)

deploydocs(
    repo = "github.com/NonlinearOscillations/HarmonicBalance.jl.git",
    push_preview = false
)
