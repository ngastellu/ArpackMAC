
module constructHtb_relaxed
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./CoordsIO.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem

nstruc = ARGS[1]
structype = ARGS[2]
posfile = expanduser("~/Desktop/simulation_outputs/MAC_structures/Ata_structures/$structype/relaxed/$(structype)n$(nstruc)_relaxed.xsf")
#strucindex = parse(Int,split(split(split(posfile,'/')[end],'-')[2], '_')[1])
println("Reading coords from file: $posfile...")
fullpos, _ = read_xsf(posfile; read_forces=false)
# fullpos = get_frame(posfile, 1)

println(size(fullpos))
const rCC::Float64 = 1.8 #max nearest neighbour distance in angstrom

py"""import numpy as np
from qcnico.remove_dangling_carbons import remove_dangling_carbons
nn = $nstruc
rCC = $rCC
pos = remove_dangling_carbons($(PyObject(fullpos)),$rCC)
"""

pos = PyArray(py"pos"o)

println("Constructing hamiltonian...")
H, ii, jj, hvals = lindbergHtb_sparse(pos,rCC, return_data=true)
println("Done!")

py"""import numpy as np
nn= $nstruc
structype = $structype
hdir = f'/Users/nico/Desktop/simulation_outputs/percolation/Ata_structures/{structype}/Hao_ARPACK/local/'
np.save(hdir + f"hvals/hvals-{nn}.npy",$(PyObject(hvals)))
np.save(hdir + f"inds/ii-{nn}.npy",$(PyObject(ii)))
np.save(hdir + f"inds/jj-{nn}.npy",$(PyObject(jj)))
"""
end
