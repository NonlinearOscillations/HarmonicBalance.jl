using FileIO


function export_results_dict(res::Result)
	# takes an instance of the Result Structure and ommits the NaN's
	# by creating a dictionary at each elelement of Result.solutions
	# which only saves the non-Nan indices and values

    dimensions_parameters = size(res.solutions)
    dimension_solutions = length(res.solutions[1])
        
    # THis is a janky way of asking whether we sweep through 1 or 2 parameters
    if length(dimensions_parameters) == 1
        M=Array{Dict}(undef, dimensions_parameters[1])
    elseif length(dimensions_parameters) == 2
        M=Array{Dict}(undef, dimensions_parameters[1],dimensions_parameters[2])
    end

    # Iterate over each of the array elements
    for iter in eachindex(res.solutions)
        A = res.solutions[iter]
        value_list= []
        key_list = []
        # ask wether the solution is NaN
        for i in 1:dimension_solutions
            if isnan(A[i][1]) == false
                push!(key_list, i)
                push!(value_list, A[i])
            end
        end
        # save a dictionary at each element of the M matrix
        M[iter] = Dict(key_list .=> value_list)
    end
    return M  
end

# Example for some result obtained by result = get_steady_states(harmonic_eq2, varied, fixed)

MM = export_results_dict(result);

save("dict_version.jld2", Dict("solution_matrix" => MM, "classes" => result.classes,  
        "swept_parameters" => result.swept_parameters, "jacobian"=>result.jacobian, 
         "num_solutions" => length(result.solutions[1,1])))



