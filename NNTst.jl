module NNTst
    include("./TightBinding.jl")
    include("./CoordsIO.jl")
    
    using .TightBinding, .CoordsIO
    using LinearAlgebra, PyCall

    fn = "data/bigMAC_10x10-64_relaxed.xsf"
    pos,_ = read_xsf(fn; read_forces=false)
    pos = pos[:,1:2]
    rCC = 1.8
    N = size(pos,1)

    py"""import numpy as np 
    inn = np.load("data/nn_inds.npy")
    dnn = np.load("data/nn_dists.npy")"""

    innpy = PyArray(py"inn"o) .+ 1
    dnnpy = PyArray(py"dnn"o)

    #println(innpy[1:5,:])
    #println(innpy[end-5:end,:])

    dnn, ii, jj = nn_pairdists_vec(pos, rCC)
    println(typeof(ii))
    inn = cat(ii,jj,dims=2)
    sinn = Set(eachrow(inn))
    sinnpy = Set(eachrow(innpy))
    sample = rand(sinn,10)
    sample2 = rand(sinnpy,10)

    for (s,s2) in zip(sample,sample2)
        println("$s  $s2")
    end

    Δ = symdiff(sinn,sinnpy)
    # println("Δ = ", Δ)
    println("inn == innpy: ", issetequal(sinn,sinnpy))

    # for ij ∈ Δ
    #     println("Inds = $ij ; dnn$ij = $(dnn[ij]) ; dnnpy$ij = $(dnnpy[ij])")
    # end

    @time nn_pairdists_vec(pos,rCC)




    # println(all([inn[n,1] <= inn[n+1,1] for n=1:(size(inn,1)-1)]))

    # idiffs = inn .!= innpy
    # wherediffs = findall(idiffs) # 1D vector of CartesianIndex objects
    # tmp_i = sortperm(wherediffs,by=x->x[1]) # sorting tuples sorts them according to their 1st element
    # wherediffs = wherediffs[tmp_i]
    # println(sum([dij[2] == 1 for dij ∈ wherediffs]))
    # println(sum([dij[2] == 2 for dij ∈ wherediffs]))

    # ndiffs = size(wherediffs,1)
    # bothdiffs = findall([wherediffs[i][1] == wherediffs[i+1][1] for i=1:ndiffs-1])
    # #println(inn[wherediffs[bothdiffs]])
    # #println(innpy[wherediffs[bothdiffs]])

    # wherediffs = map(x->CartesianIndex(x[1],x[2]+1),wherediffs)
    # println(inn[wherediffs[bothdiffs]])
    # println(innpy[wherediffs[bothdiffs]])

    # println(maximum(inn))
    # println(maximum(innpy))

    # iM = argmax(inn)[1]
    # println(inn[iM,:])

end