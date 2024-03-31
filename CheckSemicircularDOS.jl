module CheckSemicircularDOS

include("./TightBinding.jl")

using LinearAlgebra, PyCall, Statistics


nn = ARGS[1]
N = parse(Int, ARGS[2])

println("Running on structure $nn")

μ = -2.39415246
σ = 0.13658531


py"""import numpy as np
n = $nn
ii = np.load(f'hamiltonians/inds/ii-{n}.npy')
jj = np.load(f'hamiltonians/inds/jj-{n}.npy')
"""

ii = PyArray(py"ii"o)
jj = PyArray(py"jj"o)

println("Number of atoms = $N\nConstructing hamiltonian...")

nb_elems = size(ii,1)

hvals = μ .+ (σ .* randn(nb_elems))

Hgaussian = zeros((N,N))

for k=1:nb_elems
   i = ii[k] 
   j = jj[k]
   h = hvals[k]
   Hgaussian[i,j] = h
end

println("Done! Diagonalising...")

eigvals = LAPACK.syev!('N','U', Hgaussian)

println("Done!")

py"""
n = $nn
np.save(f'gaussian_hvals/gvals-{n}.npy', $(PyObject(Hgaussian)))
np.save(f'gaussian_eigvals/eigvals-{n}.npy', $(PyObject(eigvals)))
"""

end
