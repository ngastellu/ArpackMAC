module run_QuickArpackBigMAC
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./CoordsIO.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem, MPI


MPI.Init()
comm = MPI.COMM_WORLD

nprocs = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)


#posfile = expanduser(ARGS[1])
#strucindex = parse(Int,split(split(split(posfile,'/')[end],'-')[2], '_')[1])

temp = parse(Int, ARGS[1])
frame0_index = parse(Int, ARGS[2])
nframes = parse(Int, ARGS[3])
framestep = parse(Int, ARGS[4])

const trajfile::String = "../../lammps_MD/$(temp)K_no-rotate/dump_traj.xsf"
const rCC::Float64 = 1.8 #max nearest neighbour distance in angstrom


println("[$rank] Reading coords from file: $trajfile...")
flush(stdout)

nframes_proc = Int(nframes/nprocs)
frame0_proc = frame0_index + nframes_proc * framestep * rank


for n=0:nframes_proc-1

    frame_index = frame0_proc + n*framestep
    println("[$rank] ******** Working on frame $frame_index ********")
    flush(stdout)


    fullpos = get_frame(trajfile, frame_index)

    println("[$rank] $(size(fullpos))")
    flush(stdout)

    py"""import numpy as np
    from qcnico.remove_dangling_carbons import remove_dangling_carbons
    rCC = $rCC
    pos = remove_dangling_carbons($(PyObject(fullpos)),$rCC)
    """

    pos = PyArray(py"pos"o)

    println("[$rank] Constructing hamiltonian...")
    flush(stdout)

    H = lindbergHtb_sparse(pos,rCC)
    println("[$rank] Done!")
    flush(stdout)

    py"""import numpy as np
    frame = $frame_index
    np.save(f"H-{frame}.npy",$(PyObject(H)))
    """

    N = size(H,1)

    if N % 2 != 0
        error("Number of atoms needs to be even! We have N = $N.") 
    end

    nhalf = Int(N/2)
    #eps_QCFFPI = 2.7e-7
    eps_tb = 1e-7
    
    print("[$rank] Estimating eHOMO...")
    flush(stdout)
    approx_eHOMO, Rspectrum = estimate_eHOMO(H,eps_tb*100) #Rspectrum = spectral range
    print(" [$rank] Done! ")
    flush(stdout)
    println("[$rank] Estimated eHOMO = $(approx_eHOMO) eV")
    flush(stdout)
    println("[$rank] Spectral range = $Rspectrum eV")
    flush(stdout)

    print("Running one-shot Lanczos... ")
    flush(stdout)
    ε, ψ, iLUMO = LUMO_arpack_MAC(H, approx_eHOMO, Rspectrum)
    nconv = size(ε,1)
    println("[$rank] Done! Obtained $nconv eigenvalues. iLUMO = $iLUMO")
    flush(stdout)

    py"""import numpy as np
    nn = $frame_index
    ii = $iLUMO
    np.save(f"eARPACK_bigMAC_iLUMO={ii}-{nn}.npy",$(PyObject(ε)))
    np.save(f"MOs_ARPACK_bigMAC_iLUMO={ii}-{nn}.npy",$(PyObject(ψ)))    
    """
end

end