module GetExtremalMOs 

	include("./QuickArpackBigMAC.jl")
    
    using .QuickArpackBigMAC,  PyCall, Base.Filesystem, SparseArrays

    
    
    ii = parse(Int, ARGS[1])


    py"""
    import numpy as np
    from qcnico.qcffpi_io import read_MO_file, read_Hao
    from os import path

    qcffpi_dir = path.expanduser('~/scratch/bigMAC_qcffpi/')

    ii=$ii

    with open(qcffpi_dir + 'successful_runs.txt') as fo:
        for i in range(ii):
            line = fo.readline() 
    n = int(line.strip())
    pos, _ = read_MO_file(qcffpi_dir + f'sample-{n}/MO_coefs.dat')
    N = pos.shape[0]
    Hao = read_Hao(qcffpi_dir + f'sample-{n}/Hao.dat', N)
    """

    n = PyAny(py"n"o)
    pos = PyArray(py"pos"o)
    Hao = PyArray(py"Hao"o)
    
    println("**** $n ****")
    
    H = sparse(Hao)

    
    N = size(H,1)
    nhalf = Int(floor(N/2))
    #eps_QCFFPI = 2.7e-7
    eps_tb = 1e-7
    
    print("Estimating eHOMO...")
    approx_eHOMO, Rspectrum = estimate_eHOMO(H,eps_tb*100) #Rspectrum = spectral range
    print("Done! ")
    println("Estimated eHOMO = $(approx_eHOMO) eV")

    print("Running one-shot Lanczos... (occ) ")
    εocc, ψocc = kBT_arpack_MAC(H, approx_eHOMO;MO_type="occupied")
    nconv = size(ε,1)
    println("Done! Obtained $nconv eigenvalues.")

    print("Running one-shot Lanczos... (vir) ")
    εvir, ψvir = kBT_arpack_MAC(H, approx_eHOMO;MO_type="virtual")
    nconv = size(ε,1)
    println("Done! Obtained $nconv eigenvalues.")
    
    py"""
    nn = $n
    np.save(f'energies/occupied/eARPACK_occ_bigMAC-{nn}.npy', $(PyObject(εocc)))
    np.save(f'MOs/occupied/MOs_ARPACK_occ_bigMAC-{nn}.npy', $(PyObject(ψocc)))
    np.save(f'energies/virtual/eARPACK_vir_bigMAC-{nn}.npy', $(PyObject(εvir)))
    np.save(f'MOs/virtual/MOs_ARPACK_vir_bigMAC-{nn}.npy', $(PyObject(ψvir)))
    """
    


end
