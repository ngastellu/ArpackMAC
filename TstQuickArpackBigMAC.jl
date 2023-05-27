module TstQuickArpackBigMAC
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")

using .QuickArpackBigMAC, .SpectralLanczos
using Plots, LinearAlgebra, SparseArrays, PyCall

#structure_index = 64

print("Loading TB data from Python...")
@pyinclude("./macTB.py")
H = sparse(PyArray(py"tbsys.H"o))
#ε_py = PyArray(py"tbsys.energies"o)
#ψ_py = PyArray(py"tbsys.M"o)
println("Done!")

N = size(H,1)
nhalf = Int(N/2)
eps_QCFFPI = 2.7e-7
 
print("Estimating eLUMO...")
approx_eLUMO = estimate_eLUMO(H,eps_QCFFPI*100)
print("Done! ")
println("Estimated eLUMO = $(approx_eLUMO) eV")

print("Running one-shot Lanczos... ")
ε, ψ = one_shot_arpack_MAC(H, approx_eLUMO, 300.0, eps_QCFFPI)
nconv = size(ε,1)
println("Done! Obtained $nconv eigenvalues.")
resids = get_resids(ε,ψ,H)
p = plot(ε,resids,seriestype=:scatter,legend=false)
xlabel!("Energy [/t]")
ylabel!("Residuals")
display(p)
readline()




# p = plot(ε_py[nhalf+1:nhalf+nconv],ε[nhalf+1:nhalf+nconv])
# xlabel!("Exact energies / t")
# ylabel!("Lanczos energies / t")
# display(p)

end