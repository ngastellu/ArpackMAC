module QuickArpackBigMAC
include("./Gershgorin.jl")
include("./SpectralLanczos.jl")

using LinearAlgebra, SparseArrays, PyCall, Arpack, Glob
using .Gershgorin, .SpectralLanczos

export estimate_eHOMO, kBT_arpack_MAC, LUMO_arpack_MAC, rerun_LUMO_arpack_MAC, check_δN, extremal_MOs, extra_extremal_MOs


function estimate_eHOMO(H,eps_loose)
    N = size(H,1)

    # * Step 1: Roughly estimate where HOMO energy is
    λmin, λmax = spectral_bounds(H)
    println("λmin = $λmin")
    println("λmax = $λmax")
    eL_guess = 0.5 * (λmin + λmax) # initial guess is the middle of the spectrum
    nvals, eL_guess = count_evals(H, eL_guess,-1e-9) #use negative eps_shift to avoid missing HOMO

    # Define some useful constants
    Rspectrum = λmax - λmin
    ugs = Rspectrum/N #unit guess shift
    n0 = 1 + Int(floor(N/2)) #number of levels STRICTLY under HOMO 
    δN = nvals - n0
    
    # This approach might be too coarse... might never converge
    while abs(δN) > 0.1 * N #if our initial guess of HOMO energy is way off
        eL_guess = eL_guess - ugs * δN
        nvals, eL_guess = count_evals(H,eL_guess,-1e-9)
        δN = nvals - n0
    end

    # Here, we run an under-converged Lanczos to get a finer estimate of eHOMO.
    # This will tell us what sigma to choose when running the actual Lanczos
    while δN != 0
        println("δN1 = $δN")
        if δN > 0 #overerestimated eHOMO
            elanczos, nconv, _, _, _ = eigs(H,nev=δN,which=:SR,sigma=eL_guess-1e-9,ritzvec=false,tol=eps_loose,check=1)
            sort!(elanczos)
            eL_guess = elanczos[1]
        
        elseif δN < 0 #understimated eHOMO
            elanczos, nconv, _, _, _ = eigs(H,nev=-δN,which=:LR,sigma=eL_guess,ritzvec=false,tol=eps_loose,check=1)
            sort!(elanczos)
            eL_guess = elanczos[end]
        end
        nvals, eL_guess = count_evals(H,eL_guess,-1e-9)
        δN = nvals - n0
        println("δN2 = $δN")
    end
    println("δN final = $δN")
    return eL_guess, Rspectrum
end

function check_δN(H,ε;type="LUMO")
    N = size(H,1)
    if type == "LUMO"
        ntarget = Int(floor(N/2)) #nb of MOs below LUMO
    else
        ntarget = Int(floor(N/2)) -1 #nb of MOs below HOMO
    end
    nvals = size(ε,1)
    n = 0
    δN = 1
    while (δN != 0) && (n < nvals)
        if type == "LUMO"
            n+=1 
            println("ε[$n] = $(ε[n])")
            δN = count_evals(sparse(H),ε[n])[1] - ntarget
            println("n = $n\t δN = $δN")
            if δN > 0
                # δN should be 1 bc if not, it means lanczos missed an eigenpair which lies in the range 
                # of the eigenpairs we have already obtained
                println("δN = $δN > 0; exiting now n =$n (<----- should be 1)")  
                break
            end
        else
            m = nvals - n 
            println("ε[$m] = $(ε[m])")
            δN = count_evals(sparse(H),ε[m])[1] - ntarget
            println("m = $m\t δN = $δN")
            if δN > 0
                # δN should be nvals bc if not, it means lanczos missed an eigenpair which lies in the range 
                # of the eigenpairs we have already obtained
                println("δN = $δN < 0; exiting now m =$m (<----- should be $nvals)")  
                break
            end
            n+=1
        end
    end
    return δN, n
end

function kBT_arpack_MAC(H,eHOMO_guess,T=400.0,eps_lanczos=1e-8;MO_type="virtual")
    N = size(H,1)
    kB = 8.617e-5 #eV/K
    if MO_type != "virtual" && MO_type != "occupied"
        println("[kBT_arpack_MAC] Invalid value for argument MO_type: $MO_type.\nMust be 'occupied' or 'virtual'(default).\nSetting MO_type = 'virtual'.")
        MO_type = "virtual"
    end
    if MO_type == "virtual"
        eboundary = eHOMO_guess + 3.0 * kB * T #this is the maximum energy of the states we want
        nvals, _ = count_evals(H,eboundary)
        nev_req = nvals - Int(floor(N/2)) 
        println("Number of requested eigenvalues =  $(nev_req)")
        ε,ψ, _, _, _, _ = eigs(H,nev=nev_req,sigma=eHOMO_guess-1e-8,which=:LR,maxiter=10000,tol=eps_lanczos)
    else
        eboundary = eHOMO_guess + 1e-9 - 1.5 * kB * T #this is the minimum energy of the states we want; add 1e-8 to eHOMO to avoid missing iti
        println(eboundary)
        nvals, _ = count_evals(H,eboundary)
        nev_req = Int(floor(N/2)) - nvals
        println("Number of requested eigenvalues =  $(nev_req)")
        ε,ψ, _, _, _, _ = eigs(H,nev=nev_req,sigma=eHOMO_guess+1e-9,which=:SR,maxiter=10000,tol=eps_lanczos)
    end
    return ε,ψ
end

function LUMO_arpack_MAC(H, eHOMO_guess, Rspectrum; eps_lanczos=1e-9, nvals=4) #Rspectrum = spectral range of H; obtained in `estimate_eHOMO` above (from Gershgorin circle thm)
    N = size(H,1)
    ugs = Rspectrum/(2*N) #This unit step value is smaller than in `estimate_eHOMO` bc we expect |δN| < 1 (the estimate of HOMO should be quite already; no need to go too crazy with the reshifts)
    ε,ψ, _, _, _, _ = eigs(H,nev=nvals,sigma=eHOMO_guess-2 *ugs,which=:LM,maxiter=10000,tol=eps_lanczos) #take look for eigenvalues near energies decently lower than the eHOMO estimate to avoid missing HOMO and to get HOMO
    δN, iLUMO = check_δN(H,ε)
    ntries = 1
    while δN != 0 
        println("**** try nb. $ntries; δN = $δN **** ")
        if δN > 0
            ε0 = ε[1] 
        else
            ε0 = ε[end] 
        end
        if ntries % 10 == 0
            eHOMO_guess = ε0 - ugs * δN * 20 #if algo is stagnating, give it a big kick every 10 steps
        else
            eHOMO_guess = ε0 - ugs * δN
        end
        ε,ψ, _, _, _, _ = eigs(H,nev=nvals,sigma=eHOMO_guess-ugs,which=:LM,maxiter=10000,tol=eps_lanczos)
        δN, iLUMO = check_δN(H,ε)
        ntries += 1
    end
    return ε, ψ, iLUMO
end

function get_closest_eLUMO(nframe)

    existing_efiles = glob("eARPACK_bigMAC_iLUMO*npy")
    inds = [parse(Int, split(split(f, '.')[1], '-')[2]) for f in existing_efiles]
    closest_frame_file_ind = argmin(abs.(nframe .- inds))
    efile = existing_efiles[closest_frame_file_ind]
    iLUMO_last = parse(Int, split(split(efile, '=')[2], '-')[1])
    
    py"""
    import numpy as np
    energies = np.load($efile)
    """

    last_eLUMO = PyArray(py"energies"o)[iLUMO_last]

    return last_eLUMO, iLUMO_last

end

function rerun_LUMO_arpack_MAC(H, eHOMO_guess, Rspectrum, T, nframe; eps_lanczos=1e-9,nvals=4)
    N = size(H,1)
    ugs = Rspectrum/N #This unit step value is smaller than in `estimate_eHOMO` bc we expect |δN| < 1 (the estimate of HOMO should be quite already; no need to go too crazy with the reshifts)
    ε,ψ, _, _, _, _ = eigs(H,nev=nvals,sigma=eHOMO_guess-2 *ugs,which=:LM,maxiter=10000,tol=eps_lanczos) #take look for eigenvalues near energies decently lower than the eHOMO estimate to avoid missing HOMO and to get HOMO
    δN, iLUMO = check_δN(H,ε)
    ntries = 1
    py"""import numpy as np
    T = $T
    distrib_params = np.load(f'../energy_distribution_params/{T}K_edistribution_params.npy')"""
    μ_e, σ_e, emode = PyArray(py"distrib_params"o)
    d_e = Normal(μ_e, σ_e)
    d_ugs = Uniform(0.05,2.0)
    while δN != 0 && ntries < 200
        println("[frame = $nframe] **** try nb. $ntries; δN = $δN **** ")
        if δN > 0
            ε0 = ε[1] 
        else
            ε0 = ε[end] 
        end
        if ntries == 10
            eLUMO_guess, iLUMO_last = get_closest_eLUMO(nframe)
            println("[frame = $nframe]: Closest completed frame = $iLUMO_last")
        elseif ntries == 20
           eLUMO_guess = μ_e
        elseif ntries == 40
            eLUMO_guess = emode
        elseif ntries % 10 == 0 #every 10 tries, draw new HOMO guess from
            eLUMO_guess = rand(d_e)
        else
            eLUMO_guess = ε0 - ugs * δN * rand(d_ugs) #randomise shift to avoid ping-ponging between the same values
        end
        ε,ψ, _, _, _, _ = eigs(H,nev=nvals,sigma=eLUMO_guess-(ugs/10),which=:LM,maxiter=10000,tol=eps_lanczos)
        δN, iLUMO = check_δN(H,ε)
        ntries += 1
    end
    if δN != 0
        println("LUMO not converged after $ntries tries; δN = $δN. Returning all 1s.")
        ε = ones(4)
        ψ = ones(4,4)
        iLUMO=5
    end
    return ε, ψ, iLUMO
end

function extremal_MOs(H;nvals_lo=50,nvals_hi=50,eps_lanczos=1e-9)
        εlo,ψlo, _, _, _, _ = eigs(H,nev=nvals_lo,which=:SR,maxiter=100000,tol=eps_lanczos)
        εhi,ψhi, _, _, _, _ = eigs(H,nev=nvals_hi,which=:LR,maxiter=100000,tol=eps_lanczos)
        return εlo,ψlo,εhi,ψhi 
end

function extra_extremal_MOs(H, highest_elo, lowest_ehi; nvals=50,eps_lanczos=1e-9,eps_shift=1e-6)
        εlo,ψlo, _, _, _, _ = eigs(H,nev=nvals,sigma=highest_elo+eps_shift,which=:LR,maxiter=100000,tol=eps_lanczos)
        εhi,ψhi, _, _, _, _ = eigs(H,nev=nvals,sigma=lowest_ehi-eps_shift,which=:SR,maxiter=100000,tol=eps_lanczos)
        return εlo,ψlo,εhi,ψhi 
end
end