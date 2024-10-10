module GetExtremalMOs 

	include("./QuickArpackBigMAC.jl")
    
    using .QuickArpackBigMAC,  PyCall, Base.Filesystem, SparseArrays


    function get_nvals(H::SparseMatrixCSC,T::Number,motype::String)
        # Estimates the number of eigenvalues within 3*kB*T of the highest (if `motype` = 'hi) or lowest eigenvalue 
        # (if `motype` = 'lo')  of H
        kB = 8.617e-5 #eV/K
        if motype == "lo"
            e0, φ = eig(H,nev=1,which=:SR)
            eboundary = e0 + 3*kB*T
            nvals, eboundary = count_evals(H,eboundary)
        elseif motype == "hi"
            emax, φ = eig(H,nev=1,which=:LR)
            eboundary = emax - 3*kB*T
            nvals_under, eboundary = count_evals(H,eboundary)
            
            N = size(H,1)
            nvals = N - nvals_under # n eigvals are < eboundary ⟹ (N-n) eigvals ≥ eboundary
        else
            println("[get_nvals ERROR] MO type $motype is invalid! Must be either \"hi\" or \"lo\". Returning 0.")
            nvals = 0
        end
        return nvals
    end


    py"""
    import numpy as np
    evenN_lbls = np.load('evenN_lbls.npy')
    """

    lbls = PyArray(py"evenN_lbls"o)
    rCC = 1.8
    T = 400 #K

    ii = parse(Int, ARGS[1])
    n = lbls[ii]
  
    
    println("**** $n ****")
    

    natoms_file = "natoms_even.npy"
    natoms = PyArray(py"np.load($natoms_file)"o)

    ii_file = "hamiltonians/inds/ii-$n.npy"
    iii = PyArray(py"np.load($ii_file)"o)

    jj_file = "hamiltonians/inds/jj-$n.npy"
    jj = PyArray(py"np.load($jj_file)"o)

    hvals_file = "hamiltonians/hvals/hvals-$n.npy"
    hvals = PyArray(py"np.load($hvals_file)"o)

    N = natoms[ii]
 
    println("Constructing hamiltonian...")
    H = sparse(iii,jj,hvals,N,N)
    H += H'
    println("Done!")

    nvals_lo = get_nvals(H,T,"lo")
    nvals_hi = get_nvals(H,T,"hi")

    εlo,ψlo,εhi,ψhi = extremal_MOs(H;nvals_lo=nvals_lo,nvals_hi=nvals_hi)

    py"""
    nn = $n
    np.save(f'energies/lo/eARPACK_lo_bigMAC-{nn}.npy', $(PyObject(εlo)))
    np.save(f'MOs/lo/MOs_ARPACK_lo_bigMAC-{nn}.npy', $(PyObject(ψlo)))
    np.save(f'energies/hi/eARPACK_hi_bigMAC-{nn}.npy', $(PyObject(εhi)))
    np.save(f'MOs/hi/MOs_ARPACK_hi_bigMAC-{nn}.npy', $(PyObject(ψhi)))
    """
    


end
