module FFTWExt

using DSP: DSP
using FFTW: fft, fftfreq, fftshift
using Peaks: Peaks

"""
Fourier transform the timeseries of a simulation in the rotating frame and calculate the quadratures and frequencies in the non-rotating frame.
"""
function FFT(soln_u, soln_t; window=DSP.Windows.hanning)
    "Input: solution object of DifferentialEquation (positions array and corresponding time)
    Output: Fourier transform and frequencies, where window function window was used"
    w = window(length(soln_t))
    dt = soln_t[2] - soln_t[1]

    soln_tuples = Tuple.(zip(soln_u, soln_t))

    fft_u =
        length(soln_t) / sum(w) *
        [fftshift(fft(w .* [u[j] for (u, t) in soln_tuples])) for j in 1:length(soln_u[1])]
    fft_f = fftshift(fftfreq(length(soln_t), 1 / dt))

    # normalize fft_u
    return (fft_u / length(fft_f), 2 * pi * fft_f)
end

# function HarmonicBalance.FFT(soln::OrdinaryDiffEqTsit5.ODESolution; window=DSP.Windows.hanning)
#     return HarmonicBalance.FFT(soln.u, soln.t; window=window)
# end

# function FFT_analyze(fft_u::Vector{ComplexF64}, fft_f)
#     "finds peaks in the spectrum and returns corresponding frequency, amplitude and phase.
#     Frequency and phase are corrected according to Huang Dishan, Mechanical Systems and Signal Processing (1995) 9(2), 113–118
#     This correction works for a rectangular window."

#     # retaining more sigdigits gives more ''spurious'' peaks
#     max_indices, mxval = Peaks.peakproms(round.(abs.(fft_u), sigdigits=3); minprom=1)
#     Del = fft_f[2] - fft_f[1] # frequency spacing
#     A1 = abs.(fft_u)[max_indices]
#     df = zeros(length(max_indices))

#     # correction to the amplitude and phase of the peak
#     for i in 1:length(max_indices)
#         if abs.(fft_u)[max_indices[i] - 1] < abs.(fft_u)[max_indices[i] + 1]
#             A2 = abs.(fft_u)[max_indices[i] + 1]
#             df[i] = -Del / (A1[i] / A2 + 1)
#         else
#             A2 = abs.(fft_u)[max_indices[i] - 1]
#             df[i] = Del / (A1[i] / A2 + 1)
#         end
#     end
#     return 2 * pi * (fft_f[max_indices] - df),
#     A1 .* 2,
#     angle.(fft_u)[max_indices] + pi * df / Del
# end

function u_of_t(omegas_peak, As_peak, phis_peak, t)
    "Calculate us or vs as a function of time from the Fourier components."
    N = length(omegas_peak)
    u = zeros(length(t))
    for m in (Int(N / 2 - mod(N / 2, 1)) + 1):N
        u .+= As_peak[m] * cos.(omegas_peak[m] .* t .+ phis_peak[m])
    end
    return u
end

function uv_nonrotating_frame(
    omega_rot, omega_peak, A_u_peak, phi_u_peak, A_v_peak, phi_v_peak
)
    "calculates amplitudes and frequencies of the position in the nonrotating frame from the
    amplitudes and frequencies in the rotating frame."
    omega_nr = [omega_rot - omega_peak, omega_rot + omega_peak]
    u_nr =
        [
            -A_u_peak * cos(phi_u_peak) + A_v_peak * sin(phi_v_peak)
            -A_u_peak * cos(phi_u_peak) - A_v_peak * sin(phi_v_peak)
        ] ./ 2
    v_nr =
        [
            A_v_peak * cos(phi_v_peak) + A_u_peak * sin(phi_u_peak)
            A_v_peak * cos(phi_v_peak) - A_u_peak * sin(phi_u_peak)
        ] ./ 2
    return omega_nr, u_nr, v_nr
end
end # module
