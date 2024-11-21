"""
$(TYPEDEF)

Holds the three parameters of a Lorentzian peak, defined as A / sqrt((ω-ω0)² + Γ²).

# Fields
$(TYPEDFIELDS)

"""
struct Lorentzian{T<:Real}
    ω0::T
    Γ::T
    A::T
    Lorentzian(; ω0, Γ) = new{eltype(ω0)}(ω0, Γ, one(ω0)) # default peak height is 1
    function Lorentzian(ω0, Γ, A)
        type = promote_type(typeof.((ω0, Γ, A))...)
        return new{type}(convert.(type, (ω0, Γ, A))...)
    end
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
mutable struct JacobianSpectrum{T<:Real}
    peaks::Vector{Lorentzian{T}}
end
JacobianSpectrum{T}() where {T<:Real} = JacobianSpectrum{T}(Lorentzian{T}[])

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
end
