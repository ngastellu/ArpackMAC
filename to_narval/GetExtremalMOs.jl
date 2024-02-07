module GetExtremalMOs 

	include("./QuickArpackBigMAC.jl")
    
    using .QuickArpackBigMAC,  PyCall, Base.Filesystem, SparseArrays

    py"""
    import numpy as np
    evenN_lbls = np.load('evenN_lbls.npy')
    """

    lbls = PyArray(py"evenN_lbls"o)
    rCC = 1.8

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

    εlo,ψlo,εhi,ψhi = extremal_MOs(H)

    py"""
    nn = $n
    np.save(f'energies/lo/eARPACK_lo_bigMAC-{nn}.npy', $(PyObject(εlo)))
    np.save(f'MOs/lo/MOs_ARPACK_lo_bigMAC-{nn}.npy', $(PyObject(ψlo)))
    np.save(f'energies/hi/eARPACK_hi_bigMAC-{nn}.npy', $(PyObject(εhi)))
    np.save(f'MOs/hi/MOs_ARPACK_hi_bigMAC-{nn}.npy', $(PyObject(ψhi)))
    """
    


end
