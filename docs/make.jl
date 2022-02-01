push!(LOAD_PATH, "../src/")

using Documenter
using HarmonicBalance

makedocs(
	sitename="HarmonicBalance.jl",
	#format = Documenter.HTML(assets = [])
	pages = [
		"Examples" => Any[
			"examples/ex1.md"
			"examples/ex2.md"
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
