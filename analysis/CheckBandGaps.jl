module CheckBandGaps

    include("./SpectralLanczos.jl")
    include("./TightBinding.jl")
    include("./CoordsIO.jl")
    
    using .SpectralLanczos, .TightBinding, .CoordsIO, PyCall, Base.Filesystem


    structype = basename(pwd())

    if structype == "40x40"
        lbls = 1:300
    elseif structype == "tempdot6"
        lbls = 0:217
    elseif structype == "tempdot5"
        lbls = 0:216
    else
        throw(ArgumentError(structype, "not valid sAMC ensemble.")
    end

    for (k, n) ∈ enumerate(lbls)
        println("**** $n ****")
        eocc_file = expanduser("~/Desktop/simulation_outputs/percolation/40x40/eARPACK/occupied/eARPACK_bigMAC-$n.npy")
        eoccs = sort(PyArray(py"np.load($eocc_file)"o))

        evirt_file = expanduser("~/Desktop/simulation_outputs/percolation/40x40/eARPACK/virtual/eARPACK_bigMAC-$n.npy")
        evirts = sort(PyArray(py"np.load($evirt_file)"o))

        posfile = expanduser("~/Desktop/simulation_outputs/percolation/40x40/structures/bigMAC-$(n)_relaxed.xsf") 
        pos, _ = read_xsf(posfile;read_forces=false)
        py"""
        from qcnico.remove_dangling_carbons import remove_dangling_carbons
        rCC = $rCC
        pos = remove_dangling_carbons($(PyObject(pos)),rCC)
        """
        pos = PyArray(py"pos"o)
        
        println("Constructing hamiltonian...")
        H = lindbergHtb_sparse(pos,rCC)
        println("Done!")

        N = size(H,1)
        nhalf = Int(floor(N/2))

        eHOMO = eoccs[end]
        eLUMO = evirts[1]

        Δn[k,1], _ = count_evals(H, eHOMO)
        Δn[k,2], _ = count_evals(H, eLUMO)
        print('\n')

    end

    py"""
    np.save('gap_check.npy', $(PyObject(Δn)))
    """

end