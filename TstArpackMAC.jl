module TstArpackMAC
include("./ArpackMAC.jl")

using LinearAlgebra, SparseArrays, .ArpackMAC, Plots

resids, ΔE, elanczos, nmissing =  arpackMAC(64)
println("Nb of missing eigenvalues = ", nmissing)

p = plot(abs.(ΔE), resids, seriestype=:scatter,legend=false)
ylabel!("Residuals")
xlabel!("|ΔE| (TB-QCFFPI - Lanczos) [eV]")
#savefig(p,"10x10_bigMAC-64_resids_v_abs_E_err.png")
display(p)
readline()

eF = -0.131810388916
#eF = 0.0
p = plot(elanczos .- eF, resids,seriestype=:scatter,legend=false)
ylabel!("Residuals")
xlabel!("Elanczos - E_LUMO[eV]")
#savefig(p,"10x10_bigMAC-64_resids_v_E.png")
display(p)
readline()
sort!(ΔE)
println(ΔE[1])
println(ΔE[end-1])
println(ΔE[end])
end