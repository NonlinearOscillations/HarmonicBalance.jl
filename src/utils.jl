is_real(x) = abs(imag(x)) / abs(real(x)) < IM_TOL || abs(x) < 1e-70
is_real(x::Array) = is_real.(x)

flatten(a) = collect(Iterators.flatten(a))

_parse_expression(exp) = exp isa String ? Num(eval(Meta.parse(exp))) : exp
_symidx(sym::Num, args...) = findfirst(x -> isequal(x, sym), _free_symbols(args...))

tuple_to_vector(t::Tuple) = [i for i in t]
_str_to_vec(s::Vector) = s
_str_to_vec(s) = [s]

function common_type(A::Array)
    types = unique(reduce(vcat, [unique(typeof.(ps)) for ps in A]))
    return promote_type(types...)
end
function collect_parameters(sweeps::OrderedDict, fixed_parameters::OrderedDict)
    param_ranges = collect(values(sweeps))
    iter = Iterators.product(param_ranges..., values(fixed_parameters)...)
    return collect(iter)
end
function type_stable_parameters(sweeps::OrderedDict, fixed_parameters::OrderedDict)
    input_array = collect_parameters(sweeps, fixed_parameters)
    T = common_type(input_array)
    return [convert.(T, ps) for ps in input_array]
end
function parameter_type(sweeps::OrderedDict, fixed_parameters::OrderedDict)
    input_array = collect_parameters(sweeps, fixed_parameters)
    return common_type(input_array)
end

"Remove occurrences of `sweeps` elements from `fixed_parameters`."
function filter_duplicate_parameters(sweeps, fixed_parameters)
    new_params = copy(fixed_parameters)
    for par in keys(sweeps)
        delete!(new_params, par)
    end
    return new_params
end

"
Show fields of an object.
"
function show_fields(object)
    for field in fieldnames(typeof(object)) # display every field
        display(string(field))
        display(getfield(object, field))
    end
end
