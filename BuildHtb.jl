
module BuildHtb
include("./TightBinding.jl")
include("./CoordsIO.jl")

using .TightBinding, .CoordsIO
using Base.Filesystem, NPZ

structype = ARGS[1]
nx = ARGS[2]
ny = ARGS[3]

h_exists = isfile("hamiltonians/hvals/hvals-$(nx)x$(ny).npy") # if hvals exists, the Hamiltonian has already been constructed
if h_exists
	println("Hamiltonian for $(nx)x$(ny) $structype GNR already exists! My job here is done.")
else
	println("Building hamiltonian...")
    posfile = "gnr_$(structype)_$(nx)x$(ny).xyz"
    println("Reading coords from file: $posfile...")
    _, pos = PyArray(py"pos"o)
    rCC = 1.8
	H, ii, jj, hvals = lindbergHtb_sparse(pos,rCC;return_data=true)
    npzwrite("hamiltonians/hvals/hvals-$(nx)x$(ny).npy", hvals)
    npzwrite("hamiltonians/inds/ii-$(nx)x$(ny).npy", ii)
    npzwrite("hamiltonians/inds/jj-$(nx)x$(ny).npy", jj)
    println("Done!")
end

end