module NNTst
    include("./TightBinding.jl")
    include("./CoordsIO.jl")
    
    using .TightBinding, .CoordsIO
    using LinearAlgebra, PyCall, Plots

    fn = "data/bigMAC_40x40-100_relaxed.xsf"
    pos,_ = read_xsf(fn; read_forces=false)
    pos = pos[:,1:2]
    rCC = 1.8
    N = size(pos,1)

   small_pos = pos[1:100,:]

   dat = nn_pairdists(small_pos,rCC)
   dat = nn_pairdists_vec(small_pos, rCC)

#     println("***** Running nn_pairdists *****")
#    @time nn_pairdists(pos,rCC)

   println("***** Running nn_pairdists_vec *****")
   @time nn_pairdists_vec(pos,rCC)

end