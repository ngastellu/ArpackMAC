module CheckDeltaN

include("./SpectralLanczos.jl")

using SparseArrays, PyCall, .SpectralLanczos

# **** NOTA BENE: This script just checks if the LUMO obtained in `run_QuickArpackBigMAC*jl` is actually the LUMO
#                 AFTER the run is complete. It is therefore for postprocessing only. See `QuickArpackBigMAC` for
#                 the `check_δN` function that is used DURING the run. ****

frame_index = ARGS[1]

py"""import numpy as np
nn = $frame_index 
H = np.load(f'H-{nn}.npy')
ee = np.load(f'eARPACK_bigMAC-{nn}.npy')
"""

H = PyArray(py"H"o)
ee = PyArray(py"ee"o)

N = size(H)[1]
ntarget = Int(N/2)
nvals = size(ee)[1]
println("N = $N ⟹ ntarget = $ntarget")
δN = 1 #initialise δN to some random value

n = 0
while (δN != 0) && (n < nvals)
    global n+=1
    println("ee[n] = $(ee[n])")
    global δN = count_evals(sparse(H),ee[n])[1] - ntarget
    println("n = $n\t δN = $δN")
    if δN > 0
        break
    end
end

end

end