push!(LOAD_PATH, "../src/")

using Documenter
using HarmonicBalance

makedocs(
	sitename="HarmonicBalance.jl",
    modules = [HarmonicBalance],
    warnonly = true,
	format = Documenter.HTML(
		mathengine=MathJax(),
        canonical="https://nonlinearoscillations.github.io/HarmonicBalance.jl/stable/",
		assets = ["assets/favicon.ico", "assets/docs.css"]
        # size_threshold = nothing
	),
	pages = [
		"Background" => Any[
			"background/harmonic_balance.md"
			"background/stability_response.md"
		],
		"Examples" => Any[
			"examples/simple_Duffing.md"
			"examples/linear_response.md"
			"examples/time_dependent.md"
			"examples/parametron.md"
		],
		"Manual" => Any[
			"manual/entering_eom.md"
			"manual/extracting_harmonics.md"
			"manual/solving_harmonics.md"
            "manual/Krylov-Bogoliubov_method.md"
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
