module run_QuickArpackBigMAC
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .TightBinding
using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem

nstruc = ARGS[1]
structype = ARGS[2]
if structype == "40x40"
    posfile = expanduser("~/scratch/$(structype)/relaxed_structures_no_dangle/bigMAC-$(nstruc)_relaxed.xyz")
else
    posfile = expanduser("~/scratch/$(structype)/relaxed_structures_no_dangle/$(structype)n$(nstruc)_relaxed.xyz")
end
#strucindex = parse(Int,split(split(split(posfile,'/')[end],'-')[2], '_')[1])
println("Reading coords from file: $posfile...")

py"""from qcnico.coords_io import read_xyz
import numpy as np
posfile=$posfile
pos = read_xyz(posfile)
nn = $nstruc

hvals = np.load(f'hamiltonians/hvals/hvals-{nn}.npy')
ii = np.load(f'hamiltonians/inds/ii-{nn}.npy')
ii = np.load(f'hamiltonians/inds/ii-{nn}.npy')
"""

pos = PyArray(py"pos"o)
N = size(pos,1)

ii = PyArray(py"ii"o)
jj = PyArray(py"jj"o)
hvals = PyArray(py"hvals"o)

println("Constructing hamiltonian...")
H = sparse(ii, jj, hvals, N, N)
H += H'
println("Done!")

N = size(H,1)
nhalf = Int(ceil(N/2))
#eps_QCFFPI = 2.7e-7
eps_tb = 1e-7
 
print("Estimating eHOMO...")
approx_eHOMO, Rspectrum = estimate_eHOMO(H,eps_tb*100) #Rspectrum = spectral range
print("Done! ")
println("Estimated eHOMO = $(approx_eHOMO) eV")

print("Running one-shot Lanczos... ")
ε, ψ = kBT_arpack_MAC(H, approx_eHOMO;MO_type="virtual",keep_HOMO=true)
nconv = size(ε,1)
println("Done! Obtained $nconv eigenvalues.")

δn = count_evals(H, ε[1])[1] - (nhalf - 1)
println("lvldiff δn = $(δn) (should be 0)")

py"""import numpy as np
nn= $nstruc
np.save(f"/Users/nico/Desktop/simulation_outputs/percolation/40x40/eARPACK/virtual_w_HOMO/eARPACK_bigMAC-{nn}.npy",$(PyObject(ε)))
np.save(f"/Users/nico/Desktop/simulation_outputs/percolation/40x40/MOs_ARPACK/virtual_w_HOMO/MOs_ARPACK_bigMAC-{nn}.npy",$(PyObject(ψ)))
"""
end
