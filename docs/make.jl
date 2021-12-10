using Pkg
current_path = @__DIR__
Pkg.activate(current_path * "/../.");

using Documenter, HarmonicBalance

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
			]
		]
)
