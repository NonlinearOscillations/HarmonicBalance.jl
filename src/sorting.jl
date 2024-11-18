"""
    $(TYPEDSIGNATURES)

Sorts `solutions` into branches according to the method `sorting`.

`solutions` is an n-dimensional array of Vector{Vector}. Each element describes a set of solutions for a given parameter set.
The output is a similar array, with each solution set rearranged such that neighboring solution sets have the smallest Euclidean distance.

Keyword arguments
- `sorting`: the method used by `sort_solutions` to get continuous solutions branches.  The current options are `"hilbert"` (1D sorting along a Hilbert curve), `"nearest"` (nearest-neighbor sorting) and `"none"`.
- `show_progress`: Indicate whether a progress bar should be displayed.

"""
function sort_solutions(
    solutions::Solutions(T); sorting="nearest", show_progress=true
) where {T}
    sorting_schemes = ["none", "hilbert", "nearest"]
    sorting ∈ sorting_schemes ||
        error("Only the following sorting options are allowed:  ", sorting_schemes)
    sorting == "none" && return solutions
    l = length(size(solutions))
    l == 1 && return sort_1D(solutions; show_progress=show_progress)
    l == 2 && return sort_2D(solutions; sorting=sorting, show_progress=show_progress)
    return error("do not know how to solve solution which are not 1D or 2D")
end

function sort_solutions!(res::Result; sorting="nearest", show_progress=true)
    return res.solutions .= sort_solutions(
        res.solutions; sorting=sorting, show_progress=show_progress
    )
end

#####
# SOLUTION SORTING
###

"""
Removes rows and columns with given indices from the Matrix M.
The row/column indices are defined with respect to the original M!
"""
function remove_rows_columns(M::Matrix, rows::Vector{Int}, cols::Vector{Int})
    a, b = size(M)
    return M[[!in(i, rows) for i in 1:a], [!in(j, cols) for j in 1:b]]
end

remove_rows_columns(M, row::Int, col::Int) = remove_rows_columns(M, [row], [col])
function remove_rows_columns(M, pairs::Vector{CartesianIndex{2}})
    return remove_rows_columns(M, getindex.(pairs, 1), getindex.(pairs, 2))
end
remove_rows_columns(M, pair::CartesianIndex{2}) = remove_rows_columns(M, pair[1], pair[2])

"Return true if array contains any element more than once"
is_repetitive(array) = length(unique(array)) !== length(array)

"Given two sets of [u1, u2, ... un] and [v1, v2, ... vn], return a matrix of norms
such that M_ij = norm(ui - vj). NaNs (nonexistent solutions) are replaced by Inf"
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

"""
Match a to_sort vector of solutions to a set of reference vectors of solutions.
Returns a list of Tuples of the form (1, i1), (2, i2), ... such that
reference[1] and to_sort[i1] belong to the same branch
"""
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

"""
Go through a vector of solution and sort each according to Euclidean norm.
"""
function sort_1D(solns::Solutions(T); show_progress=true) where {T<:Number}
    sorted_solns = similar(solns) # preallocate
    sorted_solns[1] = sort(solns[1]; by=x -> abs.(imag(x))) # prefer real solution at first position

    if show_progress
        bar = Progress(
            length(solns); dt=1, desc="Ordering solutions into branches ...", barlen=50
        )
    end
    for i in eachindex(solns[1:(end - 1)])
        show_progress ? ProgressMeter.next!(bar) : nothing
        matched_indices = align_pair(sorted_solns[i], solns[i + 1]) # pairs of matching indices
        next_indices = getindex.(matched_indices, 2) # indices of the next solution
        sorted_solns[i + 1] = (solns[i + 1])[next_indices]
    end
    return sorted_solns
end

function hilbert_indices(solns::Solutions(T)) where {T}
    """Get mapping between 2D indexes (parameter space) and a 1D Hilbert curve"""
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

function naive_indices(solns::Solutions(T)) where {T}
    idx_pairs = []
    for i in 1:size(solns, 1)
        for j in 1:size(solns, 2)
            push!(idx_pairs, [i, j])
        end
    end
    return idx_pairs
end

function get_nn_2D(idx::Vector{Int}, Nx::Int, Ny::Int)
    "returns all neighbors from a point, including diagonal ones"
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

function sort_2D(solns::Solutions(T); sorting="nearest", show_progress=true) where {T}
    """match each 2D solution with all its surrounding neighbors, including the diagonal ones"""
    # determine a trajectory in 2D space where nodes will be visited
    if sorting == "hilbert" # propagating matching of solutions along a hilbert_curve in 2D
        idx_pairs = hilbert_indices(solns)
    elseif sorting == "nearest" # propagate matching of solutions along rows
        idx_pairs = naive_indices(solns)
    end

    # infinite solutions are ignored by the align_pair function. This trick allows a consistent ordering "propagation"
    sorted_solns = Inf .* copy(solns)
    sorted_solns[1, 1] = sort(solns[1, 1]; by=x -> abs.(imag(x))) # prefer real solution at first position

    if show_progress
        bar = Progress(
            length(idx_pairs); dt=1, desc="Ordering solutions into branches ...", barlen=50
        )
    end
    for i in 1:(length(idx_pairs) - 1)
        show_progress ? ProgressMeter.next!(bar) : nothing
        neighbors = get_nn_2D(idx_pairs[i + 1], size(solns, 1), size(solns, 2))
        reference = [sorted_solns[ind...] for ind in neighbors]
        matched_indices = align_pair(reference, solns[idx_pairs[i + 1]...]) # pairs of matching indices
        next_indices = getindex.(matched_indices, 2) # indices of the next solution
        sorted_solns[idx_pairs[i + 1]...] = (solns[idx_pairs[i + 1]...])[next_indices]
    end
    return sorted_solns
end

"Find a branch order according `classification`. Place branches where true occurs earlier first."
function find_branch_order(classification::Vector{BitVector})
    branches = [getindex.(classification, k) for k in 1:length(classification[1])] # array of branches
    indices = replace(findfirst.(branches), nothing => Inf)
    negative = findall(x -> x == Inf, indices) # branches not true anywhere - leave out
    return order = setdiff(sortperm(indices), negative)
end

find_branch_order(classification) = collect(1:length(classification[1])) # no ordering for >1D

"Order the solution branches in `res` such that close classified positively by `classes` are first."
function order_branches!(res::Result, classes::Vector{String})
    for class in classes
        order_branches!(res, find_branch_order(res.classes[class]))
    end
end

order_branches!(res::Result, class::String) = order_branches!(res, [class])

"Reorder the solutions in `res` to match the index permutation `order`."
function order_branches!(res::Result, order::Vector{Int})
    res.solutions .= _reorder_nested(res.solutions, order)
    for key in keys(res.classes)
        res.classes[key] = _reorder_nested(res.classes[key], order)
    end
end
