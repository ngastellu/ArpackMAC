module CheckBandGaps

    include("./SpectralLanczos.jl")
    include("./TightBinding.jl")
    include("./CoordsIO.jl")
    
    using .SpectralLanczos, .TightBinding, .CoordsIO, PyCall, Base.Filesystem

    py"""
    import numpy as np
    vir_lbls = np.load('vir_lbls.npy')
    """

    lbls = PyArray(py"vir_lbls"o)[1:2]
    rCC = 1.8

    ii = parse(Int, ARGS[1])
    n = lbls[ii]
    Δn = zeros(Int, 2)
    
    println("**** $n ****")
    eocc_file = "energies/occupied/eARPACK_bigMAC-$n.npy"
    eoccs = sort(PyArray(py"np.load($eocc_file)"o))

    evirt_file = "energies/virtual/eARPACK_bigMAC-$n.npy"
    evirts = sort(PyArray(py"np.load($evirt_file)"o))

    posfile = expanduser("~/scratch/clean_bigMAC/40x40/relax/no_PBC/relaxed_structures/bigMAC-$(n)_relaxed.xsf") 
    pos, _ = read_xsf(posfile;read_forces=false)
    py"""
    from qcnico.remove_dangling_carbons import remove_dangling_carbons
    rCC = $rCC
    pos = remove_dangling_carbons($(PyObject(pos)),rCC)
    """
    pos = PyArray(py"pos"o)
    
    println("Constructing hamiltonian...")
    H, ii, jj, hvals = lindbergHtb_sparse(pos,rCC;return_data=true)
    println("Done!")

    N = size(H,1)
    nhalf = Int(floor(N/2))

    eHOMO = eoccs[end]
    eLUMO = evirts[1]

    Δn[k,1], _ = count_evals(H, eHOMO)
    Δn[k,2], _ = count_evals(H, eLUMO)


    py"""
    nn = $n
    np.save(f'gap_check/dN-{nn}.npy', $(PyObject(Δn)))
    np.save(f'hamiltonians/inds/ii-{nn}.npy', $(PyObject(ii)))
    np.save(f'hamiltonians/inds/jj-{nn}.npy', $(PyObject(jj)))
    np.save(f'hamiltonians/hvals/hvals-{nn}.npy', $(PyObject(hvals)))
    """

end