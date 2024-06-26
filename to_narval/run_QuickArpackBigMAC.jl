module run_QuickArpackBigMAC
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./CoordsIO.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem

nstruc = ARGS[1]
posfile = expanduser("~/scratch/clean_bigMAC/40x40/relax/no_PBC/relaxed_structures/bigMAC-$(nstruc)_relaxed.xsf")
#strucindex = parse(Int,split(split(split(posfile,'/')[end],'-')[2], '_')[1])
println("Reading coords from file: $posfile...")
fullpos, _ = read_xsf(posfile; read_forces=false)
# fullpos = get_frame(posfile, 1)

println(size(fullpos))
const rCC::Float64 = 1.8 #max nearest neighbour distance in angstrom

py"""import numpy as np
from remove_dangling_carbons import remove_dangling_carbons
rCC = $rCC
pos = remove_dangling_carbons($(PyObject(fullpos)),$rCC)
"""

pos = PyArray(py"pos"o)

println("Constructing hamiltonian...")
H = lindbergHtb_sparse(pos,rCC)
println("Done!")

N = size(H,1)
nhalf = Int(floor(N/2))
#eps_QCFFPI = 2.7e-7
eps_tb = 1e-7
 
print("Estimating eHOMO...")
approx_eHOMO, Rspectrum = estimate_eHOMO(H,eps_tb*100) #Rspectrum = spectral range
print("Done! ")
println("Estimated eHOMO = $(approx_eHOMO) eV")

print("Running one-shot Lanczos... ")
ε, ψ = kBT_arpack_MAC(H, approx_eHOMO;MO_type="occupied")
nconv = size(ε,1)
println("Done! Obtained $nconv eigenvalues.")

py"""import numpy as np
nn= $nstruc
np.save(f"energies/occupied/eARPACK_bigMAC-{nn}.npy",$(PyObject(ε)))
np.save(f"MOs/occupied/MOs_ARPACK_bigMAC-{nn}.npy",$(PyObject(ψ)))
"""
end
