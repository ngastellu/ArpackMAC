
module run_QuickArpackBigMAC
    include("./QuickArpackBigMAC.jl")
    include("./SpectralLanczos.jl")
    include("./CoordsIO.jl")
    include("./TightBinding.jl")
    
    using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
    using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem
    
    #posfile = expanduser(ARGS[1])
    #strucindex = parse(Int,split(split(split(posfile,'/')[end],'-')[2], '_')[1])
    posfile = "/Users/nico/Desktop/simulation_outputs/percolation/40x40/structures/bigMAC-210_relaxed.xsf"
    println("Reading coords from file: $posfile...")
    fullpos, _ = read_xsf(posfile; read_forces=false)
    
    const rCC::Float64 = 1.8 #max nearest neighbour distance in angstrom
    
    py"""import numpy as np
    from coords_io import read_xyz
    from remove_dangling_carbons import remove_dangling_carbons
    rCC = $rCC
    fullpos = $fullpos 
    pos = remove_dangling_carbons(fullpos,$rCC)
    np.save("pos.npy", pos)
    """
    
    pos = PyArray(py"pos"o)
    
    println("Constructing hamiltonian...")
    H = lindbergHtb_sparse(pos,rCC)
    println("Done!")

    py"""
    np.save("H.npy", $(PyObject(H)))
    """
    
    N = size(H,1)
    nhalf = Int(floor(N/2))
    #eps_QCFFPI = 2.7e-7
    const eps_tb = 1e-9
    const T = 400
     
    print("Estimating eHOMO...")
    approx_eHOMO, Rspectrum = estimate_eHOMO(H,eps_tb*100) #Rspectrum = spectral range
    print("Done! ")
    println("Estimated eHOMO = $(approx_eHOMO) eV")
    
    print("Running one-shot Lanczos... ")
    ε, ψ = kBT_arpack_MAC(H, approx_eHOMO, T, eps_tb; MO_type="occupied")
    nconv = size(ε,1)
    δN, iHOMO = check_δN(H,ε;type="HOMO")
    println("Done! Obtained $nconv eigenvalues. iHOMO = $iHOMO")
    
    py"""import numpy as np
    ii = $iHOMO
    np.save(f"eARPACK_bigMAC_iHOMO={ii}.npy",$(PyObject(ε)))
    np.save(f"MOs_ARPACK_bigMAC_iHOMO={ii}.npy",$(PyObject(ψ)))
    """
    end
    