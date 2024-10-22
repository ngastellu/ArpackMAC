module CheckLUMO

# This script checks if the lowest eigenvalues obtained by kBT_arpack_MAC are LUMO, HOMO, or neither


include("./SpectralLanczos.jl")
include("./TightBinding.jl")

using .SpectralLanczos, PyCall, .TightBinding, SparseArrays

ensemble = ARGS[1]

if ensemble == "40x40"
    lbls=range(1,300)
elseif ensemble == "tempdot6"
    lbls=range(0,131)
elseif ensemble == "tempdot5"
    lbls=range(0,116)
end

ncheck = 4 # number of eigenvalues to check
δn = ones(Int,size(lbls,1),ncheck+1) .* -1

for (k,n) ∈ enumerate(lbls)
    print("$n: ")
    py"""
import numpy as np
n = $n
try:
    ii = np.load(f'hamiltonians/inds/ii-{n}.npy')
    jj = np.load(f'hamiltonians/inds/jj-{n}.npy')
    hvals = np.load(f'hamiltonians/hvals/hvals-{n}.npy')
    M = np.load(f'MOs/virtual/MOs_ARPACK_bigMAC-{n}.npy')
    ee = np.load(f'energies/virtual/eARPACK_bigMAC-{n}.npy')
    N = M.shape[0]
except FileNotFoundError as e:
    print(e)
    N = -1
    """
    N = convert(Int, py"N"o)

    # handle case where no NPY file was loaded by above Python code
    if N == -1
        continue
    end

    ii = PyArray(py"ii"o)
    jj = PyArray(py"jj"o) 
    hvals = PyArray(py"hvals"o)

    ee = PyArray(py"ee"o)

    H = sparse(ii, jj, hvals, N, N)
    H += H'

    δn[k,1] = n

    for j=1:ncheck
        δn[k,j+1] = count_evals(H, ee[j])[1] - Int(ceil(N/2))
        print(" $(δn[k,j+1])")
    end
    print('\n')
end

py"""
delta_n = $({PyObject(δn)})
delta_n = delta_n[delta_n[0] >= 0]]
np.save('delta_ns.npy', delta_n)
"""

end