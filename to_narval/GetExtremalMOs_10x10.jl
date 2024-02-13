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

    εlo,ψlo,εhi,ψhi = extremal_MOs(H)

    py"""
    nn = $n
    np.save(f'energies/lo/eARPACK_lo_bigMAC-{nn}.npy', $(PyObject(εlo)))
    np.save(f'MOs/lo/MOs_ARPACK_lo_bigMAC-{nn}.npy', $(PyObject(ψlo)))
    np.save(f'energies/hi/eARPACK_hi_bigMAC-{nn}.npy', $(PyObject(εhi)))
    np.save(f'MOs/hi/MOs_ARPACK_hi_bigMAC-{nn}.npy', $(PyObject(ψhi)))
    """
    


end
