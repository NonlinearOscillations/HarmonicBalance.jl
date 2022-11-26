using HarmonicBalance

dummy_soln_vector = [[rand(ComplexF64, 10) for i in 1:10] for pt in 1:10]
dummy_soln_matrix = [[rand(ComplexF64, 10) for i in 1:10] for pt in CartesianIndices(zeros(10, 10))]

HarmonicBalance.sort_1D(dummy_soln_vector, show_progress=false)
HarmonicBalance.sort_2D(dummy_soln_matrix, show_progress=false)