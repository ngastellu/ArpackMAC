
module BuildHtb
include("./TightBinding.jl")
include("./CoordsIO.jl")

using .TightBinding, .CoordsIO
using Base.Filesystem, NPZ

structype = ARGS[1]
nx = ARGS[2]
ny = ARGS[3]

h_exists = isfile("hamiltonians/hvals/hvals-$(nx)x$(ny)_pbc.npy") # if hvals exists, the Hamiltonian has already been constructed
if h_exists
	println("Hamiltonian for $(nx)x$(ny) $structype GNR already exists! My job here is done.")
else
	println("Building hamiltonian...")
    posfile = "gnr_$(structype)_$(nx)x$(ny).xyz"
    println("Reading coords from file: $posfile...")
    pos, supercell, _ = read_xyz_supercell(posfile)
    pos = pos[1:2,:]
    supercell = supercell[2:3] # ASE formats the cell s.t. PBC are along z, our coords are defined s.t. PBC are along y
    rCC = 1.8
	H, ii, jj, hvals = lindbergHtb_sparse(pos,rCC;return_data=true,cellsize=supercell)
    npzwrite("hamiltonians/hvals/hvals-$(nx)x$(ny)_pbc.npy", hvals)
    npzwrite("hamiltonians/inds/ii-$(nx)x$(ny)_pbc.npy", ii)
    npzwrite("hamiltonians/inds/jj-$(nx)x$(ny)_pbc.npy", jj)
    println("Done!")
end

end