# Saving and loading

All of the types native to `HarmonicBalance.jl` can be saved into a `.jld2` file using `save` and loaded using `load`. Most of the saving/loading is performed using the package `JLD2.jl`, with the addition of reinstating the symbolic variables in the `HarmonicBalance` namespace (needed to parse expressions used in the plotting functions) and recompiling stored functions (needed to evaluate Jacobians). As a consequence, composite objects such as `Result` can be saved and loaded with no loss of information.

The function `export_csv` saves a .csv file which can be plot elsewhere.

```@docs; canonical=false
HarmonicBalance.save
HarmonicBalance.load
HarmonicBalance.export_csv
``` 