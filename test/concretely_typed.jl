using Pkg, HarmonicBalance
Pkg.add(; url="https://github.com/gdalle/CheckConcreteStructs.jl")

using CheckConcreteStructs

all_concrete(HarmonicBalance; verbose=true)
all_concrete(HarmonicBalance.HarmonicVariable)
all_concrete(HarmonicBalance.WarmUp)
all_concrete(HarmonicBalance.TotalDegree)
# all_concrete(HarmonicBalance.Polyhedral) # does not work on CI?
all_concrete(HarmonicBalance.Result)
all_concrete(HarmonicBalance.Problem)
all_concrete(HarmonicBalance.HarmonicEquation)
all_concrete(HarmonicBalance.DifferentialEquation)
all_concrete(HarmonicBalance.HarmonicVariable)
all_concrete(HarmonicBalance.AdiabaticSweep)
all_concrete(HarmonicBalance.HarmonicVariable)

all_concrete(HarmonicBalance.LinearResponse.Lorentzian)
all_concrete(HarmonicBalance.LinearResponse.ResponseMatrix)
all_concrete(HarmonicBalance.LinearResponse.JacobianSpectrum)
