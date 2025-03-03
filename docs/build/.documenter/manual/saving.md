
# Saving and loading {#Saving-and-loading}

All of the types native to `HarmonicBalance.jl` can be saved into a `.jld2` file using `save` and loaded using `load`. Most of the saving/loading is performed using the package `JLD2.jl`, with the addition of reinstating the symbolic variables in the `HarmonicBalance` namespace (needed to parse expressions used in the plotting functions) and recompiling stored functions (needed to evaluate Jacobians). As a consequence, composite objects such as `Result` can be saved and loaded with no loss of information.

The function `export_csv` saves a .csv file which can be plot elsewhere.
<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.save-manual-saving' href='#HarmonicBalance.save-manual-saving'><span class="jlbinding">HarmonicBalance.save</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
save(filename, object)

```


Saves `object` into `.jld2` file `filename` (the suffix is added automatically if not entered). The resulting file contains a dictionary with a single entry.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/saving.jl#L1-L7)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.load-manual-saving' href='#HarmonicBalance.load-manual-saving'><span class="jlbinding">HarmonicBalance.load</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
load(filename)

```


Loads an object from `filename`. For objects containing symbolic expressions such as `HarmonicEquation`, the symbolic variables are reinstated in the `HarmonicBalance` namespace.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/saving.jl#L22-L28)

</details>

<details class='jldocstring custom-block' open>
<summary><a id='HarmonicBalance.export_csv-manual-saving' href='#HarmonicBalance.export_csv-manual-saving'><span class="jlbinding">HarmonicBalance.export_csv</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
export_csv(filename, res, branch)

```


Saves into `filename` a specified solution `branch` of the Result `res`.


[source](https://github.com/NonlinearOscillations/HarmonicBalance.jl/blob/372cbbb0e8435a5ab0ff80b9d5ec55fed51e08fd/src/saving.jl#L77-L81)

</details>


```

```

