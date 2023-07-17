module run_QuickArpackBigMAC
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./CoordsIO.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem

posfile = expanduser(ARGS[1])
#strucindex = parse(Int,split(split(split(posfile,'/')[end],'-')[2], '_')[1])
println("Reading coords from file: $posfile...")
#fullpos, _ = read_xsf(posfile; read_forces=false)
fullpos = get_frame(posfile, 1)

println(size(fullpos))
const rCC::Float64 = 1.8 #max nearest neighbour distance in angstrom

py"""import numpy as np
from remove_dangling_carbons import remove_dangling_carbons
rCC = $rCC
pos = remove_dangling_carbons($(PyObject(fullpos)),$rCC)
np.save("pos.npy", pos)
"""

pos = PyArray(py"pos"o)

println("Constructing hamiltonian...")
H = lindbergHtb_sparse(pos,rCC)
println("Done!")

N = size(H,1)
nhalf = Int(floor(N/2))
#eps_QCFFPI = 2.7e-7
eps_tb = 1e-6
 
print("Estimating eHOMO...")
approx_eHOMO, Rspectrum = estimate_eHOMO(H,eps_tb*100) #Rspectrum = spectral range
print("Done! ")
println("Estimated eHOMO = $(approx_eHOMO) eV")

print("Running one-shot Lanczos... ")
ε, ψ, iLUMO = LUMO_arpack_MAC(H, approx_eHOMO, Rspectrum)
nconv = size(ε,1)
println("Done! Obtained $nconv eigenvalues. iLUMO = $iLUMO")

py"""import numpy as np
ii = $iLUMO
np.save(f"eARPACK_bigMAC_iLUMO={ii}.npy",$(PyObject(ε)))
np.save(f"MOs_ARPACK_bigMAC_iLUMO={ii}.npy",$(PyObject(ψ)))
"""
end