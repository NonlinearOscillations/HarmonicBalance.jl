```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: HarmonicBalance.jl
  text: Efficient Floquet expansions for nonlinear driven systems
  tagline: A Julia suite for nonlinear dynamics using harmonic balance
  actions:
    - theme: brand
      text: Getting Started
      link: /introduction
    - theme: alt
      text: Tutorials
      link: /tutorials
    - theme: alt
      text: View on GitHub
      link: https://github.com/NonlinearOscillations/HarmonicBalance.jl
  image:
    src: /logo.svg
    alt: HarmonicBalance.jl

features:
  - title: Non-equilibrium steady states
    details: Compute the all stationary states in a one or two-dimensional parameter sweep.
  - title: Linear Response
    details: Explore the fluctuations on top of the steady states.
  - title: Limit Cycles
    details: Find limit cycles involving many frequencies.
---
```

## How to Install HarmonicBalance.jl?

It is easy to install HarmonicBalance.jl as we are registered in the Julia General registry.
You can simply run the following command in the Julia REPL:
```julia
julia> using Pkg
julia> Pkg.add("HarmonicBalance")
```

If you want to use the latest unreleased version of HarmonicBalance.jl, you can run the following
command: (in most cases the released version will be same as the version on github)
```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/NonlinearOscillations/HarmonicBalance.jl")
```