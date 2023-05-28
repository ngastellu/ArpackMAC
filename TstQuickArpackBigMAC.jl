module TstQuickArpackBigMAC
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./CoordsIO.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
using Plots, LinearAlgebra, SparseArrays, PyCall

posfile = "data/bigMAC_40x40-100_relaxed.xsf"
pos, _ = read_xsf(posfile; read_forces=false)
println(size(pos))
rCC = 1.8 #max nearest neighbour distance in angstrom

H = lindbergHtb_sparse(pos,rCC)

N = size(H,1)
nhalf = Int(N/2)
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