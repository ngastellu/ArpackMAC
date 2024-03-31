module DenseTB

include("./TightBinding.jl")

using LinearAlgebra, PyCall, .TightBinding

nn = ARGS[1]
N = parse(Int, ARGS[2])

println("Running on structure $nn")
println("Loading entries and indices...")
py"""import numpy as np
n = $nn
ii = np.load(f'hamiltonians/inds/ii-{n}.npy')
jj = np.load(f'hamiltonians/inds/jj-{n}.npy')
hvals = np.load(f'hamiltonians/hvals/hvals-{n}.npy')
"""
ii = PyArray(py"ii"o)
jj = PyArray(py"jj"o)
hvals = PyArray(py"hvals"o)

println("Done! Constructing hamiltonian...")

H = sparse(ii,jj,hvals,N,N)
H += H'
H = Matrix(H)

println("Done! Diagonalising...")

eigvals = LAPACK.syev!('N','U', H)

println("Done!")

py"""
n = $nn
np.save(f'dense_tb_eigvals/eigvals-{n}.npy', $(PyObject(eigvals)))
"""

end
