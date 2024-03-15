module run_QuickArpackBigMAC
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./CoordsIO.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem

nstruc = ARGS[1]
structype = ARGS[2]
posfile = expanduser("~/scratch/$(structype)/$(structype)n$(nstruc).xyz")
#strucindex = parse(Int,split(split(split(posfile,'/')[end],'-')[2], '_')[1])
println("Reading coords from file: $posfile...")

const rCC::Float64 = 1.8 #max nearest neighbour distance in angstrom

py"""import numpy as np
from remove_dangling_carbons import remove_dangling_carbons
from qcnico.coords_io import read_xyz

nn = $nstruc
rCC = $rCC
pos = read_xyz($posfile)
pos = remove_dangling_carbons(pos,$rCC)
np.save(f'coords/coords-{nn}.npy', pos)
"""

pos = PyArray(py"pos"o)

println("Constructing hamiltonian...")
H, ii, jj, hvals = lindbergHtb_sparse(pos,rCC,return_data=true)
println("Done!")

py"""import numpy as np
nn = $nstruc
np.save(f'hamiltonians/inds/ii-{nn}.npy', $(PyObject(ii)))
np.save(f'hamiltonians/inds/jj-{nn}.npy', $(PyObject(jj)))
np.save(f'hamiltonians/hvals/hvals-{nn}.npy', $(PyObject(hvals)))
"""

N = size(H,1)
nhalf = Int(floor(N/2))
#eps_QCFFPI = 2.7e-7
eps_tb = 1e-7
 
print("Estimating eHOMO...")
approx_eHOMO, Rspectrum = estimate_eHOMO(H,eps_tb*100) #Rspectrum = spectral range
print("Done! ")
println("Estimated eHOMO = $(approx_eHOMO) eV")

print("Running one-shot Lanczos for occupied states... ")
εocc, ψocc = kBT_arpack_MAC(H, approx_eHOMO;MO_type="occupied")
nconv = size(εocc,1)
println("Done! Obtained $nconv eigenvalues.")

print("Running one-shot Lanczos for virtual states... ")
εvir, ψvir = kBT_arpack_MAC(H, approx_eHOMO;MO_type="virtual")
nconv = size(εvir,1)
println("Done! Obtained $nconv eigenvalues.")

py"""import numpy as np
nn= $nstruc
np.save(f"energies/occupied/eARPACK_bigMAC-{nn}.npy",$(PyObject(εocc)))
np.save(f"MOs/occupied/MOs_ARPACK_bigMAC-{nn}.npy",$(PyObject(ψocc)))
np.save(f"energies/virtual/eARPACK_bigMAC-{nn}.npy",$(PyObject(εvir)))
np.save(f"MOs/virtual/MOs_ARPACK_bigMAC-{nn}.npy",$(PyObject(ψvir)))
"""
end
