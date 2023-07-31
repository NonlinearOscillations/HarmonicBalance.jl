export @maybethread
export is_real
export *
export @variables
export d
export show

macro maybethread(loop)
    if Threads.nthreads()>1
      quote Threads.@threads $(Expr(loop.head,
                               Expr(loop.args[1].head, esc.(loop.args[1].args)...),
                               esc(loop.args[2]))); end
    else
      # @warn "running single threaded"
      quote $(esc(loop)); end
    end
end

import Base: show, display;

import Base: ComplexF64, Float64; export ComplexF64, Float64
ComplexF64(x::Complex{Num}) = ComplexF64(Float64(x.re) + im*Float64(x.im))
Float64(x::Complex{Num}) = Float64(ComplexF64(x))
Float64(x::Num) = Float64(x.val)

# default global settings
export IM_TOL
IM_TOL::Float64 = 1E-6
function set_imaginary_tolerance(x::Float64)
   @eval(IM_TOL::Float64 = $x)
end

is_real(x) = abs(imag(x)) / abs(real(x)) < IM_TOL::Float64 || abs(x) < 1e-70
is_real(x::Array) = is_real.(x)

    # Symbolics does not natively support complex exponentials of variables
    import Base: exp
    exp(x::Complex{Num}) = x.re.val == 0 ? exp(im*x.im.val) : exp(x.re.val + im*x.im.val)
