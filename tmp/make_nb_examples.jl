using Literate

using HarmonicBalance, Plots
default(; fmt=:png)

### Process examples
# Always rerun examples
const EXAMPLES_IN = joinpath(@__DIR__, "..", "docs", "src", "examples")
const OUTPUT_NB_DIR = @__DIR__

examples = filter!(file -> file[(end - 2):end] == ".jl", readdir(EXAMPLES_IN; join=true))

function preprocess(content)
    sub = SubstitutionString("""
                             """)
    content = replace(content, r"^# # [^\n]*"m => sub; count=1)

    # remove VSCode `##` block delimiter lines
    content = replace(content, r"^##$."ms => "")
    return content
end

for example in examples
    Literate.notebook(example, OUTPUT_NB_DIR; documenter=false, execute=true)
end

for example in examples
    Literate.notebook(example, OUTPUT_NB_DIR; documenter=false, execute=true)
end
