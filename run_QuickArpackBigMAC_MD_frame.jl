module run_QuickArpackBigMAC
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./CoordsIO.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem

#posfile = expanduser(ARGS[1])
#strucindex = parse(Int,split(split(split(posfile,'/')[end],'-')[2], '_')[1])

temp = parse(Int, ARGS[1])
frame_index = parse(Int, ARGS[2])
trajfile = "../../../lammps_MD/$(temp)K_no-rotate/dump_traj.xsf"

println("Reading coords from file: $trajfile...")
fullpos = get_frame(trajfile, frame_index)

println(size(fullpos))
const rCC::Float64 = 1.8 #max nearest neighbour distance in angstrom

py"""import numpy as np
from qcnico.remove_dangling_carbons import remove_dangling_carbons
rCC = $rCC
pos = remove_dangling_carbons($(PyObject(fullpos)),$rCC)
"""

pos = PyArray(py"pos"o)

println("Constructing hamiltonian...")
H = lindbergHtb_sparse(pos,rCC)
println("Done!")

py"""import numpy as np
frame = $frame_index
np.save(f"H-{frame}.npy",$(PyObject(H)))
"""

N = size(H,1)

if N % 2 != 0
    error("Number of atoms needs to be even! We have N = $N.") 
end

nhalf = Int(N/2)
#eps_QCFFPI = 2.7e-7
eps_tb = 1e-7
 
print("Estimating eLUMO...")
approx_eLUMO, Rspectrum = estimate_eLUMO(H,eps_tb*100) #Rspectrum = spectral range
print("Done! ")
println("Estimated eLUMO = $(approx_eLUMO) eV")
println("Spectral range = $Rspectrum eV")

print("Running one-shot Lanczos... ")
ε, ψ, iLUMO = LUMO_arpack_MAC(H, approx_eLUMO, Rspectrum)
nconv = size(ε,1)
println("Done! Obtained $nconv eigenvalues. iLUMO = $iLUMO")

py"""import numpy as np
ii = $iLUMO
np.save(f"eARPACK_bigMAC_iLUMO={ii}.npy",$(PyObject(ε)))
np.save(f"MOs_ARPACK_bigMAC_iLUMO={ii}.npy",$(PyObject(ψ)))
"""
end