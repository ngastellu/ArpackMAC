module ArpackTest
    include("./TightBinding.jl")
    using Arpack, SparseArrays, LinearAlgebra, .TightBinding

    H = sparse(constructHtb_dense(100,-1))

    println(typeof(eigs(H,sigma=1.01)))



end