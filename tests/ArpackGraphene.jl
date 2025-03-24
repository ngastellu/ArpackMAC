module ArpackGraphene
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./CoordsIO.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
using PyCall, LinearAlgebra, SparseArrays, Plots

gnr_lens = [10 13 18 19 27 33 40 48]
const nlens = length(gnr_lens)
const rCCgraphene = 1.421
const eps_tb = 1e-6
best_resids = zeros(nlens)
avg_resids = zeros(nlens)
worst_resids = zeros(nlens)

for (k,n) ∈ enumerate(gnr_lens)
    println("****** n = $n ******")
    println("Reading coords...")
    posfile = expanduser("~/Desktop/simulation_outputs/graphene/gnr_zigzag_$(n)x17.xyz")
    pos = read_xyz(posfile)
    py"""import numpy as np
         from remove_dangling_carbons import remove_dangling_carbons
         rCC = $rCC
         pos = remove_dangling_carbons($(PyObject(fullpos)),$rCC)
    """
    pos = PyArray(py"pos"o)
    println("Natoms = $(size(pos))")
    println("Constructing hamiltonian...")
    H = lindbergHtb_sparse(pos,rCCgraphene)
    print("Estimating eLUMO...")
    approx_eLUMO = estimate_eLUMO(H,eps_tb*100)
    print("Done! eLUMO = $approx_eLUMO eV")
    println("Running one-shot Lanczos...")
    ε, ψ = one_shot_arpack_MAC(H, approx_eLUMO, 300.0, eps_tb)
    resids = get_resids(ε,ψ,H)
    sort!(resids)
    best_resids[k] = resids[1]
    worst_resids[k] = resids[end]
    avg_resids[k] = sum(resids)/length(resids)
end

p = plot(gnr_lens,best_resids, seriestype=:scatter,label="Best")
plot!(gnr_lens,worst_resids, seriestype=:scatter,label="Worst")
plot!(gnr_lens,avg_resids, seriestype=:scatter,label="Mean")
xlabel!("Length of GNR")
ylabel!("Residuals")
display(p)
end