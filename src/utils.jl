is_real(x) = abs(imag(x)) / abs(real(x)) < IM_TOL::Float64 || abs(x) < 1e-70
is_real(x::Array) = is_real.(x)

flatten(a) = collect(Iterators.flatten(a))
