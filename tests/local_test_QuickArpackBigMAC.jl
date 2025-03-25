module run_QuickArpackBigMAC
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")


using .QuickArpackBigMAC, .SpectralLanczos
using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem

#posfile = expanduser(ARGS[1])
#strucindex = parse(Int,split(split(split(posfile,'/')[end],'-')[2], '_')[1])

py"""
import numpy as np
H = np.load('H-10009.npy')
"""

H = sparse(PyArray(py"H"o))

N = size(H,1)

if N % 2 != 0
    error("Number of atoms needs to be even! We have N = $N.") 
end

nhalf = Int(N/2)
#eps_QCFFPI = 2.7e-7
eps_tb = 1e-7
 
print("Estimating eHOMO...")
approx_eHOMO, Rspectrum = estimate_eHOMO(H,eps_tb*100) #Rspectrum = spectral range
print("Done! ")
println("Estimated eHOMO = $(approx_eHOMO) eV")
println("Spectral range = $Rspectrum eV")

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