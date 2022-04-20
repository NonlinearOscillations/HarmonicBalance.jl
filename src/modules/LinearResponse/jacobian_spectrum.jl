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


"Converts an 2-component vector belonging to a VDP variable into the FT of the natural variable."
 function _VDP_to_spectrum(λ, eigvec::Vector; ω::Float64)
    u,v = eigvec
    peak1 = (1 + 2*(imag(u)*real(v) - real(u)*imag(v)) ) * Lorentzian(ω0=ω - imag(λ) , Γ=real(λ))
    peak2 = (1 + 2*(real(u)*imag(v) - real(v)*imag(u)) ) * Lorentzian(ω0=ω + imag(λ) , Γ=real(λ))
    JacobianSpectrum([peak1, peak2])
end


#   Returns the spectra of all variables in `res` for `index` of `branch`.
function JacobianSpectrum(res::Result; index::Int, branch::Int)

    res.classes["stable"][index][branch] || error("\nThe solution is unstable - it has no JacobianSpectrum!\n")

    solution_dict = get_single_solution(res, branch=branch, index=index)

    VDPs = res.problem.eom.variables # fetch the vector of HarmonicVariable
    λs, vs = eigen(res.jacobian(solution_dict))

    solution_dict = get_single_solution(res, index=index, branch=branch)

    # blank JacobianSpectrum for each variable
    all_spectra = Dict{Num, JacobianSpectrum}([[nat_var, JacobianSpectrum([])] for nat_var in getfield.(VDPs, Symbol("natural_variable"))])
        
    for (j, λ) in enumerate(λs)
        v = vs[:, j]
        # we know how to convert each VDP variable into an excitation -> break down the eigenvector into its constituent VDP pairs
        for (i, VDP_var) in enumerate(VDPs)
            v_pair = v[[2*i-1, 2*i]] # fetch the relevant part of the Jacobian eigenvector
            ωrot = Float64(substitute(VDP_var.ω, solution_dict))
            # v_pair is associated to a natural variable -> this variable gets Lorentzian peaks
            this_spec =  norm(v_pair) * _VDP_to_spectrum(λ, v_pair, ω=ωrot)

            all_spectra[VDP_var.natural_variable] = add_peak(all_spectra[VDP_var.natural_variable], this_spec) 
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
    C = Array{Float64, 2}(undef,  length(Ω_range)-1, length(X)-1)

    for ij in CartesianIndices(C)
        C[ij] = abs(evaluate(spectra[ij[2]][nat_var], Ω_range[ij[1]] - offset[ij[2]]))
    end
    x_mat = x_scale * hcat([x*ones(length(Ω_range)) for x in X]...)
    y_mat = y_scale * hcat([Ω_range for j=1:length(X)]...)
    C = logscale ? log.(C) : C

    PyPlot.pcolormesh(x_mat, y_mat, C)
    xlabel(Latexify.latexify(string(first(keys(res.swept_parameters)))), fontsize=24);

    y_label = y_offset=="0.0" ? "noise " * latexify("ω") : "noise " * latexify("ω") * " - " * latexify(y_offset)
    ylabel(y_label, fontsize=24, fontname="Times");
    return C
end