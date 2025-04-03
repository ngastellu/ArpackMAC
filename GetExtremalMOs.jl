module GetExtremalMOs 

	include("./QuickArpackBigMAC.jl")
	include("./SpectralLanczos.jl")
    
    using .QuickArpackBigMAC, .SpectralLanczos, PyCall, Base.Filesystem, SparseArrays, Arpack


    function get_nvals(H::SparseMatrixCSC,T::Number,motype::String)
        # Estimates the number of eigenvalues within 3*kB*T of the highest (if `motype` = 'hi) or lowest eigenvalue 
        # (if `motype` = 'lo')  of H
        kB = 8.617e-5 #eV/K
        if motype == "lo"
            e0, φ = eigs(H,nev=1,which=:SR)
            eboundary = e0[1] + 3*kB*T
            nvals, eboundary = count_evals(H,eboundary)
        elseif motype == "hi"
            emax, φ = eigs(H,nev=1,which=:LR)
            eboundary = emax[1] - 3*kB*T
            nvals_under, eboundary = count_evals(H,eboundary)
            
            N = size(H,1)
            nvals = N - nvals_under # n eigvals are < eboundary ⟹ (N-n) eigvals ≥ eboundary
        else
            println("[get_nvals ERROR] MO type $motype is invalid! Must be either \"hi\" or \"lo\". Returning 0.")
            nvals = 0
        end
        return nvals
    end



	# ---------- MAIN ----------
	
    
    rCC = 1.8
    T = 400 #K 
	n = ARGS[1]
	structype = ARGS[2]
	tsoft = split(structype, '-')[2]
	posfile = expanduser("~/scratch/clean_bigMAC/final_sAMC_structures/$(structype)/sAMC$tsoft-$n.xyz")

    
	
	println("Parsing atomic coords...")

	py"""from qcnico.coords_io import read_xyz
	import numpy as np
	posfile=$posfile
	pos = read_xyz(posfile)
	"""	
	pos = PyArray(py"pos"o)
	N = size(pos,1)

    println("Done! N = $N")
	
	if N % 2 == 1
		println("!!! WARNING: odd nb of atoms !!!")
	end



    ii_file = "hamiltonians/inds/ii-$n.npy"
    iii = PyArray(py"np.load($ii_file)"o)

    jj_file = "hamiltonians/inds/jj-$n.npy"
    jj = PyArray(py"np.load($jj_file)"o)

    hvals_file = "hamiltonians/hvals/hvals-$n.npy"
    hvals = PyArray(py"np.load($hvals_file)"o)
 
    println("Constructing hamiltonian...")
    H = sparse(iii,jj,hvals,N,N)
    H += H'
    println("Done!")

#    nvals_lo = get_nvals(H,T,"lo")
#    println("nvals lo = $nvals_lo")
#    nvals_hi = get_nvals(H,T,"hi")
#    println("nvals hi = $nvals_hi")
	
	nvals_hi = 100
	nvals_lo = 100

    εlo,ψlo,εhi,ψhi = extremal_MOs(H;nvals_lo=nvals_lo,nvals_hi=nvals_hi)

    

    py"""nn = $n
stype = $structype
np.save(f'energies/lo/eARPACK_lo_{stype}-{nn}.npy', $(PyObject(εlo)))
np.save(f'MOs/lo/MOs_ARPACK_lo_{stype}-{nn}.npy', $(PyObject(ψlo)))
np.save(f'energies/hi/eARPACK_hi_{stype}-{nn}.npy', $(PyObject(εhi)))
np.save(f'MOs/hi/MOs_ARPACK_hi_{stype}-{nn}.npy', $(PyObject(ψhi)))"""
end
