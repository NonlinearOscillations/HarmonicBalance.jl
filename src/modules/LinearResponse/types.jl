"""
$(TYPEDEF)

Holds the three parameters of a Lorentzian peak, defined as A / sqrt((ω-ω0)² + Γ²).

# Fields
$(TYPEDFIELDS)

"""
struct Lorentzian
    ω0::Float64
    Γ::Float64
    A::Float64
    Lorentzian(;ω0::Float64, Γ::Float64) = new(ω0, Γ, 1) # default peak height is 1
    Lorentzian(ω0, Γ, A) = new(ω0, Γ, A)
 end


 """
 $(TYPEDEF)
 
 Holds a set of `Lorentzian` objects belonging to a variable.
 
 # Fields
 $(TYPEDFIELDS)

 # Constructor
 ```julia
 JacobianSpectrum(res::Result; index::Int, branch::Int)
 ``` 

 """
 mutable struct JacobianSpectrum
    peaks::Vector{Lorentzian}
 end


"""
$(TYPEDEF)

Holds the compiled response matrix of a system.

# Fields
$(TYPEDFIELDS)

"""
struct ResponseMatrix
    """The response matrix (compiled)."""
    matrix::Matrix{Function}
    """Any symbolic variables in `matrix` to be substituted at evaluation."""
    symbols::Vector{Num}
    """The frequencies of the harmonic variables underlying `matrix`. These are needed to transform
    the harmonic variables to the non-rotating frame."""
    variables::Vector{HarmonicVariable}

    ResponseMatrix(matrix, symbols, variables) = new(matrix, symbols, variables)
end