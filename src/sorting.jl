
using BijectiveHilbert
using LinearAlgebra
import Distances


"""
    $(TYPEDSIGNATURES)

Sorts `solutions` into branches according to the method `sorting`.

`solutions` is an n-dimensional array of Vector{Vector}. Each element describes a set of solutions for a given parameter set.
The output is a similar array, with each solution set rearranged such that neighboring solution sets have the smallest Euclidean distance.

"""
function sort_solutions(solutions::Array; sorting="hilbert")
    sorting_schemes = ["none", "hilbert", "nearest"]
    sorting ∈ sorting_schemes || error("Only the following sorting options are allowed:  ", sorting_schemes)
    sorting == "none" && return solutions
     l = length(size(solutions))
     l == 1 && return sort_1D(solutions)
     l == 2 && return sort_2D(solutions, sorting=sorting)
     error("do not know how to solve solution which are not 1D or 2D")
end


function sort_solutions!(solutions::Result; sorting="hilbert")
    solutions.solutions = sort_solutions(solutions.solutions, sorting=sorting)
end

#####
# SOLUTION SORTING
###

"""
Removes rows and columns with given idices from the Matrix M.
The row/column indices are defined with respect to the original M!
"""
function remove_rows_columns(M::Matrix, rows::Vector{Int64}, cols::Vector{Int64})
    a,b = size(M)
    M[[!in(i, rows) for i in 1:a], [!in(j, cols) for j in 1:b]]
end

remove_rows_columns(M, row::Int64, col::Int64) = remove_rows_columns(M, [row], [col])
remove_rows_columns(M, pairs::Vector{CartesianIndex{2}}) = remove_rows_columns(M, getindex.(pairs, 1), getindex.(pairs, 2))
remove_rows_columns(M, pair::CartesianIndex{2}) = remove_rows_columns(M, pair[1], pair[2])


"Return true if array contains any element more than once"
is_repetitive(array) = length(unique(array)) !== length(array)




"Given two sets of [u1, u2, ... un] and [v1, v2, ... vn], return a matrix of norms
such that M_ij = norm(ui - vj). NaNs (nonexistent solutions) are replaced by Inf"
function get_distance_matrix(ref::Vector{SteadyState}, to_sort::Vector{SteadyState})
    length(ref) != length(to_sort) && error("trying to align two solutions of unequal length!")
    distances = Distances.pairwise(Distances.Euclidean(), ref, to_sort)
    replace!(distances, NaN => Inf) # findall does not work with NaN but works with Inf
end


"Match each solution from to_sort to a closest partner from refs"
function get_distance_matrix(refs::Vector{Vector{SteadyState}}, to_sort::Vector{SteadyState})
    distances = [get_distance_matrix(ref, to_sort) for ref in refs]
    lowest_distances = similar(distances[1])
    for (i, el) in enumerate(lowest_distances)
        lowest_distances[i] = min([d[i] for d in distances]...)
    end
    lowest_distances
end



"""
Match a to_sort vector of solutions to a set of reference vector of solutions.
Returns a list of Tuples of the form (1, i1), (2, i2), ... such that
reference[1] and to_sort[i1] belong to the same branch
"""
function align_pair(reference::Vector{Vector{SteadyState}}, to_sort::Vector{SteadyState})
   
    distances = get_distance_matrix(reference, to_sort)
    n = length(to_sort)
    cartesians = [(j,i) for i in 1:n for j in 1:n]
    sorted_cartesians = cartesians[sortperm(vec(distances))]

    matched = falses(n)
    matched_ref = falses(n)

    sorted = Array{Tuple{Int64, Int64}, 1}(undef, n)
    
    for i in 1:length(cartesians)
        j,k = sorted_cartesians[i]
        if !matched[k] && !matched_ref[j]
            matched[k] = true
            matched_ref[j] = true
            sorted[j] = (j,k)
        end
    end

    return sorted

end


align_pair(ref::Vector{SteadyState}, to_sort::Vector{SteadyState}) = align_pair([ref], to_sort)


"""
Go through a vector of solution and sort each according to Euclidean norm.
"""
function sort_1D(solns::Vector{Vector{SteadyState}})
    sorted_solns = similar(solns) # preallocate
    sorted_solns[1] = sort(solns[1], by= x->abs.(imag(x))) # prefer real solution at first position

    bar = Progress(length(solns), 1, "Ordering solutions into branches ...")
    for (i, soln) in enumerate(solns[1:end-1])
        next!(bar)
        matched_indices = align_pair(sorted_solns[i], solns[i+1]) # pairs of matching indices
        next_indices = getindex.(matched_indices, 2) # indices of the next solution
        sorted_solns[i+1] = (solns[i+1])[next_indices]
    end
    sorted_solns
end

function hilbert_indices(solns::Matrix{Vector{Vector{ComplexF64}}})
    """Get mapping between 2D indexes (parameter space) and a 1D Hilbert curve"""
    Lx,Ly = size(solns) 
    mapping = [] # compute mapping between Hilbert indices and 2Ds
    for j in 1:Ly # length of parameter sweep 1
        for i in 1:Lx # length of parameter sweep 2
            X = [i, j]
            h = encode_hilbert(Simple2D(Int), X)
            X .= 0
            push!(mapping,(h=>decode_hilbert!(Simple2D(Int), X, h)))
        end 
    end
    idx_pairs = [el[2] for el in sort(mapping)] # sort along the Hilbert curve. Now we can iterate over these indexes
end

function naive_indices(solns::Matrix{Vector{Vector{ComplexF64}}})
    idx_pairs = []
    for i in 1:size(solns,1)
        for j in 1:size(solns,2)
            push!(idx_pairs,[i,j])
        end
    end
    idx_pairs
end

function get_nn_2D(idx::Vector{Int64},Nx::Int64,Ny::Int64)
    "returns all neighbors from a point, including diagonal ones"
    x, y = idx[1],idx[2]
    max_n = 1
    neighbors = []
    for x2 in x-max_n:x+max_n
        for y2 in y-max_n:y+max_n
            if (0<x<=Nx) && (0<y<=Ny) && (x != x2 || y != y2) && (1 <= x2 <= Nx) && (1 <= y2 <= Ny) 
                push!(neighbors,[x2,y2])
            end
        end
    end
    neighbors
end


function sort_2D(solns::Matrix{Vector{Vector{ComplexF64}}}; sorting="nearest") 
    """match each 2D solution with all its surrounding neighbors, including the diagonal ones"""
     # determine a trajectory in 2D space where nodes will be visited
    if sorting=="hilbert" # propagating matching of solutions along a hilbert_curve in 2D
        idx_pairs = hilbert_indices(solns)
    elseif sorting=="nearest" # propagate matching of solutions along rows
        idx_pairs = naive_indices(solns)
    end

    sorted_solns = Inf.*copy(solns)  # infinite solutions are ignored by the align_pair function. This trick allows a consistent ordering "propagation"
    sorted_solns[1,1] = sort(solns[1,1], by= x->abs.(imag(x))) # prefer real solution at first position

    bar = Progress(length(idx_pairs), 1, "Ordering solutions into branches ...")
    for i in 1:length(idx_pairs)-1 
        next!(bar)
        neighbors =  get_nn_2D(idx_pairs[i+1],size(solns,1),size(solns,2))
        reference = [sorted_solns[ind...] for ind in neighbors]
        matched_indices = align_pair(reference, solns[idx_pairs[i+1]...]) # pairs of matching indices
        next_indices = getindex.(matched_indices, 2) # indices of the next solution
        sorted_solns[idx_pairs[i+1]...] = (solns[idx_pairs[i+1]...])[next_indices]
    end 
    sorted_solns
end

