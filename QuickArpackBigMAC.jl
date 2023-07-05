module QuickArpackBigMAC
include("./Gershgorin.jl")
include("./SpectralLanczos.jl")

using LinearAlgebra, SparseArrays, PyCall, Arpack
using .Gershgorin, .SpectralLanczos

export estimate_eLUMO, one_shot_arpack_MAC


function estimate_eLUMO(H,eps_loose)

    N = size(H,1)

    # * Step 1: Roughly estimate where LUMO energy is
    λmin, λmax, indiv_bounds = spectral_bounds(H)
    println("λmin = $λmin")
    println("λmax = $λmax")
    eL_guess = 0.5 * (λmin + λmax) # initial guess is the middle of the spectrum
    nvals, eL_guess = count_evals(H, eL_guess,-1e-9) #use negative eps_shift to avoid missing LUMO

    # Define some useful constants
    Rspectrum = λmax - λmin
    ugs = Rspectrum/N #unit guess shift
    n0 = Int(floor(N/2)) #number of levels STRICTLY under LUMO 
    δN = nvals - n0
    
    # This approach might be too coarse... might never converge
    while abs(δN) > 0.1 * N #if our initial guess of LUMO energy is way off
        eL_guess = eL_guess + ugs * δN
        nvals, eL_guess = count_evals(H,eL_guess,-1e-9)
        δN = nvals - n0
    end

    # Here, we run an under-converged Lanczos to get a finer estimate of eLUMO.
    # This will tell us what sigma to choose when running the actual Lanczos
    while δN != 0
        if δN > 0 #overerestimated eLUMO
            elanczos, nconv, _, _, _ = eigs(H,nev=δN,which=:SR,sigma=eL_guess,ritzvec=false,tol=eps_loose,check=1)
            sort!(elanczos)
            eL_guess = elanczos[1]
        
        elseif δN < 0 #understimated eLUMO
            elanczos, nconv, _, _, _ = eigs(H,nev=-δN,which=:LR,sigma=eL_guess,ritzvec=false,tol=eps_loose,check=1)
            sort!(elanczos)
            eL_guess = elanczos[end]
        end
        nvals, eL_guess = count_evals(H,eL_guess,-1e-9)
        δN = nvals - n0
    end
    return eL_guess
end

function one_shot_arpack_MAC(H,eLUMO_guess,T=300.0,eps_lanczos=2e-7)
    N = size(H,1)
    kB = 8.617e-5 #eV/K
    emax = eLUMO_guess + 1.5 * kB * T
    nvals, _ = count_evals(H,emax)
    nev_req = nvals - Int(floor(N/2)) 
    println("Number of requested eigenvalues =  $(nev_req)")
    ε,ψ, _, _, _, _ = eigs(H,nev=nev_req,sigma=eLUMO_guess-1e-8,which=:LR,maxiter=10000,tol=eps_lanczos)
    return ε,ψ
end

end