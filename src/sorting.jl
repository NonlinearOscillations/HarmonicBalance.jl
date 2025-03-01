"""
    $(TYPEDSIGNATURES)

Sorts `solutions` into branches according to the specified `sorting` method.

`solutions` is an n-dimensional array of `Vector{Vector}`. Each element describes a set of
solutions for a given parameter set. The output is a similar array, with each solution set
rearranged such that neighboring solution sets have the smallest Euclidean distance.

The `sorting` keyword argument specifies the method used to get continuous solution branches.
Options are `"hilbert"` (1D sorting along a Hilbert curve), `"nearest"` (nearest-neighbor sorting),
and `"none"`. The `show_progress` keyword argument indicates whether a progress bar should be displayed.
"""
function sort_solutions(
    solutions::Solutions(T); sorting="nearest", show_progress=true
) where {T}
    sorting_schemes = ["none", "hilbert", "nearest"]
    sorting ∈ sorting_schemes ||
        error("Only the following sorting options are allowed:  ", sorting_schemes)
    sorting == "none" && return solutions
    l = length(size(solutions))
    l == 1 && return sort_1D(solutions; show_progress)
    l == 2 && return sort_2D(solutions; sorting, show_progress)
    return error("do not know how to solve solution which are not 1D or 2D")
end

"""
    $(TYPEDSIGNATURES)

Sorts the solutions in `res` in-place according to the specified `sorting` method.

`res` is a `Result` object containing the solutions to be sorted. The `sorting` keyword argument
specifies the method used to get continuous solution branches. Options are `"hilbert"`, `"nearest"`,
and `"none"`. The `show_progress` keyword argument indicates whether a progress bar should be displayed.
"""
function sort_solutions!(res::Result; sorting="nearest", show_progress=true)
    return res.solutions .= sort_solutions(res.solutions; sorting, show_progress)
end

#####
# SOLUTION SORTING
###

"Removes specified rows and columns from the matrix `M`."
function remove_rows_columns(M::Matrix, rows::Vector{Int}, cols::Vector{Int})
    a, b = size(M)
    return M[[!in(i, rows) for i in 1:a], [!in(j, cols) for j in 1:b]]
end

remove_rows_columns(M, row::Int, col::Int) = remove_rows_columns(M, [row], [col])
function remove_rows_columns(M, pairs::Vector{CartesianIndex{2}})
    return remove_rows_columns(M, getindex.(pairs, 1), getindex.(pairs, 2))
end
remove_rows_columns(M, pair::CartesianIndex{2}) = remove_rows_columns(M, pair[1], pair[2])

"Checks if the array contains any element more than once."
is_repetitive(array) = length(unique(array)) !== length(array)

"Given two sets of solutions, returns a matrix of Euclidean distances between them."
function get_distance_matrix(
    ref::Vector{SteadyState(T)}, to_sort::Vector{SteadyState(T)}
) where {T}
    length(ref) != length(to_sort) &&
        error("trying to align two solutions of unequal length!")
    distances = Distances.pairwise(Distances.Euclidean(), ref, to_sort)
    return replace!(distances, NaN => Inf) # findall does not work with NaN but works with Inf
end

"Match each solution from to_sort to a closest partner from refs"
function get_distance_matrix(refs::Solutions(T), to_sort::Vector{SteadyState(T)}) where {T}
    distances = map(ref -> get_distance_matrix(ref, to_sort), refs)
    lowest_distances = similar(distances[1])
    for idx in CartesianIndices(lowest_distances)
        lowest_distances[idx] = minimum(x[idx] for x in distances)
    end
    return lowest_distances
end

"Matches each solution from `to_sort` to the closest solution in `reference`."
function align_pair(reference, to_sort::Vector{SteadyState(T)}) where {T}
    distances = get_distance_matrix(reference, to_sort)
    n = length(to_sort)
    sorted_cartesians = CartesianIndices(distances)[sortperm(vec(distances))]

    matched = falses(n)
    matched_ref = falses(n)

    sorted = Vector{CartesianIndex}(undef, n)

    for idx in sorted_cartesians
        j, k = idx[1], idx[2]
        if !matched[k] && !matched_ref[j]
            matched[k] = true
            matched_ref[j] = true
            sorted[j] = idx
        end
    end

    return sorted
end

"Sorts 1D solutions according to the Euclidean norm."
function sort_1D(solns::Solutions(T); show_progress=true) where {T<:Number}
    sorted_solns = similar(solns) # preallocate
    # prefer real solution at first position
    sorted_solns[1] = sort(solns[1]; by=x -> abs.(imag(x)))

    if show_progress
        bar = Progress(
            length(solns); dt=1, desc="Ordering solutions into branches ...", barlen=50
        )
    end
    for i in eachindex(solns[1:(end - 1)])
        show_progress ? ProgressMeter.next!(bar) : nothing
        matched_indices = align_pair(sorted_solns[i], solns[i + 1])
        next_indices = getindex.(matched_indices, 2)
        sorted_solns[i + 1] = (solns[i + 1])[next_indices]
    end
    return sorted_solns
end

"""Get mapping between 2D indexes (parameter space) and a 1D Hilbert curve"""
function hilbert_indices(solns::Solutions(T)) where {T}
    Lx, Ly = size(solns)
    mapping = [] # compute mapping between Hilbert indices and 2Ds
    for j in 1:Ly # length of parameter sweep 1
        for i in 1:Lx # length of parameter sweep 2
            X = [i, j]
            h = encode_hilbert(Simple2D(Int), X)
            X .= 0
            push!(mapping, (h => decode_hilbert!(Simple2D(Int), X, h)))
        end
    end
    # sort along the Hilbert curve. Now we can iterate over these indexes
    return idx_pairs = [el[2] for el in sort(mapping)]
end

"Generates a list of index pairs for a 2D array of solutions."
function naive_indices(solns::Solutions(T)) where {T}
    idx_pairs = []
    for i in 1:size(solns, 1)
        for j in 1:size(solns, 2)
            push!(idx_pairs, [i, j])
        end
    end
    return idx_pairs
end

"Returns all neighbors of a point in a 2D grid, including diagonal ones."
function get_nn_2D(idx::Vector{Int}, Nx::Int, Ny::Int)
    x, y = idx[1], idx[2]
    max_n = 1
    neighbors = []
    for x2 in (x - max_n):(x + max_n)
        for y2 in (y - max_n):(y + max_n)
            if (0 < x <= Nx) &&
                (0 < y <= Ny) &&
                (x != x2 || y != y2) &&
                (1 <= x2 <= Nx) &&
                (1 <= y2 <= Ny)
                push!(neighbors, [x2, y2])
            end
        end
    end
    return neighbors
end

"""
    $(TYPEDSIGNATURES)

Sorts 2D solutions according to the specified `sorting` method.

`solns` is a 2D array of solutions. The `sorting` keyword argument specifies the method used
to get continuous solution branches.
Options are `"hilbert"` and `"nearest"`. The `show_progress` keyword argument indicates
whether a progress bar should be displayed.
"""
function sort_2D(solns::Solutions(T); sorting="nearest", show_progress=true) where {T}
    """match each 2D solution with all its surrounding neighbors, including the diagonal ones"""
    # determine a trajectory in 2D space where nodes will be visited
    if sorting == "hilbert" # propagating matching of solutions along a hilbert_curve in 2D
        idx_pairs = hilbert_indices(solns)
    elseif sorting == "nearest" # propagate matching of solutions along rows
        idx_pairs = naive_indices(solns)
    end

    # infinite solutions are ignored by the align_pair function.
    # This trick allows a consistent ordering "propagation"
    sorted_solns = Inf .* copy(solns)
    # prefer real solution at first position
    sorted_solns[1, 1] = sort(solns[1, 1]; by=x -> abs.(imag(x)))

    if show_progress
        bar = Progress(
            length(idx_pairs); dt=1, desc="Ordering solutions into branches ...", barlen=50
        )
    end
    for i in 1:(length(idx_pairs) - 1)
        show_progress ? ProgressMeter.next!(bar) : nothing
        neighbors = get_nn_2D(idx_pairs[i + 1], size(solns, 1), size(solns, 2))
        reference = [sorted_solns[ind...] for ind in neighbors]
        matched_indices = align_pair(reference, solns[idx_pairs[i + 1]...])
        next_indices = getindex.(matched_indices, 2)
        sorted_solns[idx_pairs[i + 1]...] = (solns[idx_pairs[i + 1]...])[next_indices]
    end
    return sorted_solns
end

"Finds the order of solution branches based on the classification."
function find_branch_order(classification::Vector{BitVector})
    branches = [getindex.(classification, k) for k in 1:length(classification[1])]
    indices = replace(findfirst.(branches), nothing => Inf)
    negative = findall(x -> x == Inf, indices) # branches not true anywhere - leave out
    return order = setdiff(sortperm(indices), negative)
end

find_branch_order(classification) = collect(1:length(classification[1]))

"""
Orders the solution branches in `res` such that those classified positively by `classes`
are first.
"""
function order_branches!(res::Result, classes::Vector{String})
    for class in classes
        order_branches!(res, find_branch_order(res.classes[class]))
    end
end
order_branches!(res::Result, class::String) = order_branches!(res, [class])
function order_branches!(res::Result, order::Vector{Int})
    res.solutions .= _reorder_nested(res.solutions, order)
    for key in keys(res.classes)
        res.classes[key] = _reorder_nested(res.classes[key], order)
    end
end
