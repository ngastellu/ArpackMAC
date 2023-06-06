module run_QuickArpackBigMAC
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./CoordsIO.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem

posfile = expanduser(ARGS[1])
strucindex = parse(Int,split(split(split(posfile,'/')[end],'-')[2], '_')[1])
println("Reading coords from file: $posfile...")
fullpos, _ = read_xsf(posfile; read_forces=false)

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
eps_tb = 1e-6
 
print("Estimating eLUMO...")
approx_eLUMO = estimate_eLUMO(H,eps_tb*100)
print("Done! ")
println("Estimated eLUMO = $(approx_eLUMO) eV")

print("Running one-shot Lanczos... ")
ε, ψ = one_shot_arpack_MAC(H, approx_eLUMO, 300.0, eps_tb)
nconv = size(ε,1)
println("Done! Obtained $nconv eigenvalues.")

py"""import numpy as np
nn = $strucindex
np.save(f"eARPACK_bigMAC-{nn}.npy",$(PyObject(ε)))
np.save(f"MOs_ARPACK_bigMAC-{nn}.npy",$(PyObject(ψ)))
"""
end