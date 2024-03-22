module CheckSemicircularDOS

include("./TightBinding.jl")

using LinearAlgebra, PyCall, statistics


nn = ARGS[1]

μ = 0
σ = 1.0


py"""import numpy as np
n = $nn
ii = np.load(f'hamiltonians/inds/ii-{n}.npy')
jj = np.load(f'hamiltonians/inds/jj-{n}.npy')
"""

ii = PyArray(py"ii"o)
jj = PyArray(py"jj"o)

# estimate N using max index of nonzero elements
Ni = maximum(ii)
Nj = maximum(jj)
N = maximum((Ni,Nj))

nb_elems = size(ii,1)

hvals = μ + σ * randn(nb_elems)

Hgaussian = zeros((N,N))

for k=1:nb_elems
   i = ii[k] 
   j = jj[k]
   h = hvals[k]
   Hgaussian[i,j] = h
end

eigvals = LAPACK.syev!('N','U', Hgaussian)

py"""
n = $nn
np.save(f'gaussian_hvals/gvals-{n}.npy', $(PyObject(Hgaussian)))
np.save(f'gaussian_eigvals/eigvals-{n}.npy', $(PyObject(eigvals)))
"""

end