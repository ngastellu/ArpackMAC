module constructHtb_unrelaxed
include("./QuickArpackBigMAC.jl")
include("./SpectralLanczos.jl")
include("./CoordsIO.jl")
include("./TightBinding.jl")

using .QuickArpackBigMAC, .SpectralLanczos, .CoordsIO, .TightBinding
using LinearAlgebra, SparseArrays, PyCall, Base.Filesystem

nstruc = ARGS[1]
structype = ARGS[2]
posfile = expanduser("~/Desktop/simulation_outputs/MAC_structures/Ata_structures/$structype/$(structype)n$(nstruc).xyz")
const rCC::Float64 = 1.8 #max nearest neighbour distance in angstrom

println("Reading coords from file: $posfile...")

py"""import numpy as np
from remove_dangling_carbons import remove_dangling_carbons
from qcnico.coords_io import read_xyz

nn = $nstruc
rCC = $rCC
pos = read_xyz($posfile)
pos = remove_dangling_carbons(pos,$rCC)
"""

pos = PyArray(py"pos"o)

println("Constructing hamiltonian...")
H, ii, jj, hvals = lindbergHtb_sparse(pos,rCC, return_data=true)
println("Done!")

py"""import numpy as np
nn= $nstruc
structype = $structype
hdir = f'/Users/nico/Desktop/simulation_outputs/percolation/Ata_structures/{structype}/Hao_ARPACK_unrelaxed/local/'
np.save(hdir + f"hvals/hvals-{nn}.npy",$(PyObject(hvals)))
np.save(hdir + f"inds/ii-{nn}.npy",$(PyObject(ii)))
np.save(hdir + f"inds/jj-{nn}.npy",$(PyObject(jj)))
"""
end
