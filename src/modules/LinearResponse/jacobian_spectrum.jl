export plot_jacobian_spectrum


# multiply a peak by a number.
 function *(number::Float64, peak::Lorentzian) # multiplication operation
    Lorentzian(peak.ω0, peak.Γ, peak.A*number)
 end


 *(number::Float64, s::JacobianSpectrum) = JacobianSpectrum([number * peak for peak in s.peaks])


function show(io::IO, s::JacobianSpectrum)
    peaks = sort(s.peaks, by=x->x.ω0)
    print(io, "Lorentzian peaks (central frequency ω0, linewidth Γ): \n")
    for peak in peaks
        @printf(io, "%.3e * L(ω0 = %.6e, Γ = %.3e)\n", peak.A, peak.ω0, peak.Γ)
    end
end


function show(io::IO, spectra::Dict{Num,JacobianSpectrum})
    for var in keys(spectra)
        print("VARIABLE: ", var, "\n")
        show(spectra[var])
        print("\n")
    end
end


 "Adds a peak p to JacobianSpectrum s."
 function add_peak(s::JacobianSpectrum, p::Lorentzian)
    JacobianSpectrum(cat(s.peaks, p, dims=1))
 end


 "Adds all peaks from s2 to s1."
function add_peak(s1::JacobianSpectrum, s2::JacobianSpectrum)
    for p in s2.peaks
        s1 = add_peak(s1, p)
    end
    s1
end


"Gives the numerical value of a peak at ω."
evaluate(peak::Lorentzian, ω::Float64) = peak.A / sqrt(( (peak.ω0 - ω)^2 + (peak.Γ)^2 ))

"Gives the numerical value of a JacobianSpectrum at ω"
evaluate(s::JacobianSpectrum, ω::Float64) = sum([evaluate(p, ω) for p in s.peaks])

evaluate(spectra::Dict{Num, JacobianSpectrum}, ω::Float64) = Dict([[var, evaluate(spectra[var], ω)] for var in keys(spectra)])


"Take a pair of harmonic variable u,v and an eigenvalue λ and eigenvector eigvec_2d of the Jacobian to generate corresponding Lorentzians.
    IMPORTANT: The eigenvetor eigen_2d contains only the two components of the full eigenvector which correspond to the u,v pair in question."
 function _pair_to_peaks(λ, eigvec_2d::Vector; ω::Float64)
    u,v = eigvec_2d
    peak1 = (1 + 2*(imag(u)*real(v) - real(u)*imag(v)) ) * Lorentzian(ω0=ω - imag(λ) , Γ=real(λ))
    peak2 = (1 + 2*(real(u)*imag(v) - real(v)*imag(u)) ) * Lorentzian(ω0=ω + imag(λ) , Γ=real(λ))
    JacobianSpectrum([peak1, peak2])
end


"Return the indices of a-type variables (zero harmonics) from hvars."
_get_as(hvars::Vector{HarmonicVariable}) = findall(x -> isequal(x.type, "a"), hvars)


#   Returns the spectra of all variables in `res` for `index` of `branch`.
function JacobianSpectrum(res::Result; index::Int, branch::Int)

    res.classes["stable"][index][branch] || error("\nThe solution is unstable - it has no JacobianSpectrum!\n")

    solution_dict = get_single_solution(res, branch=branch, index=index)

    hvars = res.problem.eom.variables # fetch the vector of HarmonicVariable
    λs, vs = eigen(res.jacobian(solution_dict))

    # convert OrderedDict to Dict - see Symbolics issue #601
    solution_dict = Dict(get_single_solution(res, index=index, branch=branch))

    # blank JacobianSpectrum for each variable
    all_spectra = Dict{Num, JacobianSpectrum}([[nvar, JacobianSpectrum([])] for nvar in getfield.(hvars, :natural_variable)])
        
    for (j, λ) in enumerate(λs)
        eigvec = vs[:, j] # the eigenvector

        # 2 peaks for each pair of uv variables
        for pair in _get_uv_pairs(hvars)
            u,v  = hvars[pair]
            eigvec_2d = eigvec[pair] # fetch the relevant part of the Jacobian eigenvector
            ωnum = Float64(substitute(u.ω, solution_dict)) # the harmonic (numerical now) associated to this harmonic variable

            # eigvec_2d is associated to a natural variable -> this variable gets Lorentzian peaks
            peaks =  norm(eigvec_2d) * _pair_to_peaks(λ, eigvec_2d, ω=ωnum)

            all_spectra[u.natural_variable] = add_peak(all_spectra[u.natural_variable], peaks) 
        end

        # 1 peak for a-type variable
        for a_idx in _get_as(hvars)
            a = hvars[a_idx]
            eigvec_1d = eigvec[a_idx]
            peak = 2 * norm(eigvec_1d) * Lorentzian(ω0= abs(imag(λ)) , Γ=real(λ))
            all_spectra[a.natural_variable] = add_peak(all_spectra[a.natural_variable], peak) 
        end
    end
    #_simplify_spectra!(all_spectra) # condense similar peaks
    all_spectra
end


"Condense peaks with similar centers and linewidths in JacobianSpectrum s."
function _simplify_spectrum(s::JacobianSpectrum)
    ω0s = getfield.(s.peaks, Symbol("ω0"))
    new_spectrum = JacobianSpectrum([])
    for ω in unique(ω0s)
        peaks = filter(x -> x.ω0 == ω, s.peaks) # all peaks centered at ω
        length(unique(getfield.(peaks, Symbol("Γ")))) > 1 && return(s)  # there are peaks of multiple linewidths at the same frequency => not simplifiable
        total_height = sqrt(sum(getfield.(peaks, Symbol("A")) .^2 )) # new height is the sum of squares
        new_spectrum = add_peak(new_spectrum, total_height * Lorentzian(ω0=ω,Γ=peaks[1].Γ) )
    end
    new_spectrum
end


"Simplify every JacobianSpectrum in a dictionary assigning variables to spectra."
function _simplify_spectra!(spectra::Dict{Num, JacobianSpectrum})
    for var in keys(spectra)
        spectra[var] = _simplify_spectrum(spectra[var])
    end
end 


"""
    plot_jacobian_spectrum(res::Result,
     nat_var::Num;
      Ω_range,
       branch::Int,
        y_offset::String="0.0",
          x_scale=1.0,
           y_scale=1.0,
            logscale=false)

Make a 1D plot of the response spectrum of `res` for the natural variable `nat_var`.
Performs one matrix diagonalization for each element of `Ω_range`.
This method is faster than `plot_response` but results in errors where the noise frequency
is far from the frequency of the harmonic variables.
"""
function plot_jacobian_spectrum(res::Result, nat_var::Num; Ω_range, branch::Int, y_offset::String="0.0", x_scale=1.0, y_scale=1.0, logscale=false)
    

    length(size(res.solutions)) != 1 && error("1D plots of not-1D datasets are usually a bad idea.")
    stability = classify_branch(res, branch, "stable") # boolean array
    !any(stability) && error("Cannot generate a spectrum - no stable solutions!")

    X = Vector{Float64}(collect(values(res.swept_parameters))[1][stability])

    offset = Vector{Float64}(getindex.(transform_solutions(res, y_offset), branch))[stability] 

    # only get spectra of the stable points!
    spectra = [JacobianSpectrum(res, branch=branch, index = i) for i in (1:length(res.solutions))[stability]]
    C = Array{Float64, 2}(undef,  length(Ω_range), length(X))

    bar = Progress(length(CartesianIndices(C)), 1, "Diagonalizing the Jacobian for each solution ... ", 50)
    for ij in CartesianIndices(C)
        C[ij] = abs(evaluate(spectra[ij[2]][nat_var], Ω_range[ij[1]] - offset[ij[2]]))
        next!(bar)
    end

    x_mat = x_scale * hcat([x*ones(length(Ω_range)) for x in X]...)
    y_mat = y_scale * hcat([Ω_range for j=1:length(X)]...)
    C = logscale ? log.(C) : C

    PyPlot.pcolormesh(x_mat, y_mat, C)
    xlabel(Latexify.latexify(string(first(keys(res.swept_parameters)))), fontsize=24);

    y_label = y_offset=="0.0" ? "noise " * latexify("ω") : "noise " * latexify("ω") * " - " * latexify(y_offset)
    ylabel(y_label, fontsize=24, fontname="Times");
    return X, Ω_range, C
end