is_real(x) = abs(imag(x)) / abs(real(x)) < IM_TOL || abs(x) < 1e-70
is_real(x::Array) = is_real.(x)

flatten(a) = collect(Iterators.flatten(a))

_parse_expression(exp) = exp isa String ? Num(eval(Meta.parse(exp))) : exp
_symidx(sym::Num, args...) = findfirst(x -> isequal(x, sym), _free_symbols(args...))

tuple_to_vector(t::Tuple) = [i for i in t]
_str_to_vec(s::Vector) = s
_str_to_vec(s) = [s]
