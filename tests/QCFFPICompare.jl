module QCFFPICompare
include("./SpectralLanczos.jl")
using PyCall, SparseArrays, Arpack, .SpectralLanczos

export compare_arpack_vs_qcffpi, read_energies, ekBT_ind



function get_Hao(Natoms::Int,structure_index::Int)
py"""from qcnico.qcffpi_io import read_Hao
N = $Natoms
n = $structure_index
#H = read_Hao('Hao_gnr_zigzag_5x200_TB_rescaled_unit.dat',N)
H = read_Hao(f'Hao_bigMAC_10x10-{n}.dat', N)
"""

H = PyArray(py"H"o)
end

function read_energies(structure_index)
    Ha2eV = 27.2114
    fn = "/Users/nico/Desktop/simulation_outputs/qcffpi_data/orbital_energies/orb_energy_bigMAC_10x10-$(structure_index).dat"
    #fn = "/Users/nico/Desktop/simulation_outputs/qcffpi_data/PPPvTB/orbital_energies/energies_gnr_zigzag_5x200_TB_rescaled_unit.dat"
    map(x -> parse(Float64,split(x)[2])*Ha2eV, readlines(fn))
end

function ekBT_ind(energies,T=300.0,lumo_ref=false)
    kB = 8.6217e-5
    kBT = kB * T
    N = size(energies,1)
    if lumo_ref
        elumo = energies[1+Int(floor(N/2))]
        return elumo, maximum(findall((energies .- elumo) .<= 1.5*kBT))
    else
        eF = 0.5 * (energies[Int(floor(N/2))] + energies[1+Int(floor(N/2))])
        return eF, maximum(findall((energies .- eF) .<= 1.5*kBT))
    end
end

function compare_arpack_vs_qcffpi(k,structure_index=64)
    eps_qcff_Ha = 1e-8
    Ha2eV = 27.2114
    eps_qcff_eV = eps_qcff_Ha * Ha2eV 

    eTB = read_energies(structure_index)
    N = size(eTB,1)
    H = sparse(get_Hao(N,structure_index))
    eF, ikt = ekBT_ind(eTB,300,true)
    println("ikt = $ikt, e[ikt] = $(eTB[ikt])")
    println("eF = $eF")
    nvals = ikt - Int(floor(N/2))
    elanczos, ϕ, _, niter, _, _ = eigs(H,nev=nvals,sigma=eF+3e-8,which=:LR,tol=eps_qcff_eV,check=1)
    sort_eigenpairs!(elanczos,ϕ,100000) #bs last argument
    eTB = eTB[1+Int(floor(N/2)):ikt]
    
    # Check how far-off the Lanczos energies are from their 'real values'
    ΔE = eTB .- elanczos

    # Get residuals of Ritz vectors to see how well-converged everything is
    resids = get_resids(elanczos,ϕ,H)

    # Count the number of missing eigenvalues in the spectral range obtained by our Lanczos run
    emax = elanczos[end]
    nM, M = count_evals(H, emax)
    emin = elanczos[1]
    nm, m = count_evals(H,emin,-1e-9) #use a negative shift to widen the sampled range
    if m != emin
        println("Shifted emin down by $(emin-m) eV when counting eigenvalues.")
    end
    if M != emax
        println("Shifted emax up by $(M-emax) eV when counting eigenvalues.")
    end
    nvals_expected = nM - nm
    nvals_missing = nvals_expected - size(elanczos,1)


    return resids, ΔE, elanczos, nvals_missing
end


end