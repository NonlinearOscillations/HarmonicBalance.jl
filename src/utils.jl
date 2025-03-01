is_real(x) = abs(imag(x)) / abs(real(x)) < IM_TOL || abs(x) < 1e-70
is_real(x::Array) = is_real.(x)
is_zero(x::Vector{ComplexF64}; atol=1e-6) = all(isapprox.(x, complex(0.0); atol))

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

function promote_types(sweeps::OrderedDict, fixed_parameters::OrderedDict{K}) where {K}
    param_ranges = collect(values(sweeps))
    iter = Iterators.product(param_ranges..., values(fixed_parameters)...)
    types = Iterators.flatten(unique(unique(typeof.(ps)) for ps in iter))
    T = promote_type(types...)

    promoted_sweeps = OrderedDict{K,Vector{T}}((keys(sweeps) .=> values(sweeps))...)
    promoted_fixed = OrderedDict{K,T}(fixed_parameters)
    return promoted_sweeps, promoted_fixed
end

"Remove occurrences of `sweeps` elements from `fixed_parameters`."
function filter_duplicate_parameters(sweeps, fixed_parameters)
    new_params = copy(fixed_parameters)
    for par in keys(sweeps)
        delete!(new_params, par)
    end
    return new_params
end

"Show fields of an object."
function show_fields(object)
    for field in fieldnames(typeof(object)) # display every field
        display(string(field))
        display(getfield(object, field))
    end
end

""" Project the array `a` into the real axis, warning if its contents are complex. """
function _realify(a::Array{T}; warning="") where {T<:Number}
    warned = false
    a_real = similar(a, typeof(real(a[1])))
    for i in eachindex(a)
        if !isnan(a[i]) && !warned && !is_real(a[i])
            @warn "Values with non-negligible complex parts have
            been projected on the real axis! " * warning
            warned = true
        end
        a_real[i] = real(a[i])
    end
    return a_real
end
_realify(a::Array{Real}; warning="") = a
