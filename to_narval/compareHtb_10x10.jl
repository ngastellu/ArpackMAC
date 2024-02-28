module compareHtb 

	include("./QuickArpackBigMAC.jl")
    include("./TightBinding.jl")
    
    using .QuickArpackBigMAC, .TightBinding,   PyCall, Base.Filesystem
    
    ii = parse(Int, ARGS[1])


    py"""
    import numpy as np
    from qcnico.qcffpi_io import read_MO_file, read_Hao
    from os import path

    qcffpi_dir = path.expanduser('~/scratch/bigMAC_qcffpi/')

    ii=$ii
    rCC = 1.8

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
    
    Hqcffpi = sparse(Hao)
    Harpack, ii, jj, hvals = lindbergHtb_sparse(pos, rCC; return_data=true)

    tol = 1e-13
    ΔH = droptol!(abs.(Hqcffpi - Harpack), tol)

    py"""
    nn = $n
    np.save(f'hamiltonians/inds/ii-{nn}.npy', $(PyObject(ii)))
    np.save(f'hamiltonians/inds/jj-{nn}.npy', $(PyObject(jj)))
    np.save(f'hamiltonians/hvals/hvals-{nn}.npy', $(PyObject(hvals)))
    np.save(f'diff_hamiltonians/diff-{nn}.npy', $(PyObject(ΔH)))
    """  
end
