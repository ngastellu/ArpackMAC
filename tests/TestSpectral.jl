module TestSpectral
    include("./SpectralLanczos.jl")
    include("./TightBinding.jl")

    using .SpectralLanczos, .TightBinding
    using SparseArrays, LinearAlgebra, Plots

    N = 1000
    t = -1
    A = sparse(constructHtb_dense(N, t))

    #for i=1:N
    #    A[i,i] = 0.5
    #end

    #println(A)
    println(count_evals(A,0))
    alpha = 0.0
    beta = 2.0

    E = (-2*t) .* cos.((2*pi/N).*(1:N))
    dE = abs.(E[2:end] - E[1:end-1])
    non_degen_diff = minimum(dE)

    evals, evecs, alpha, beta, last_mu = spectral_shift_lanczos(A,alpha,beta,50,1e-9,1e-12,non_degen_diff)
    print("\n\n----------------------- Post-processing -----------------------\n\n")
    nC = size(evals,1)
    nC2 = size(evecs,2)
    println("nC == nC2: ", nC == nC2)
    nA, _ = count_evals(A,alpha)
    nmu, _ = count_evals(A,last_mu)
    nC_expected = nmu - nA
    println("Expected $nC_expected Ritz pairs; got $nC ====> $(nC - nC_expected) duplicates.")
    dotprods = abs.(evecs' * evecs)
    dotprods[diagind(dotprods)] .= 0
    maxdp = maximum(dotprods)
    println(maxdp)
    dup_inds = map( Tuple,findall(dotprods .== maxdp))
    println(dup_inds)
    for tup in dup_inds
        println(typeof(tup))
        i, j = tup
        println("λ[$(i)] = $(evals[i]); λ[$(j)] = $(evals[j])")
    end

    dotprods = zeros(Int(nC * (nC-1)/2))
    Δλ = zeros(Int(nC * (nC-1)/2))
    n = 0
    for i=1:nC
        for j=1:i
            if i == j
               continue
            else
                global n += 1
                #println(n)
                Δλ[n] = abs(evals[i] - evals[j])
                dotprods[n] = abs(dot(evecs[:,i],evecs[:,j]))
            end
        end
    end
    display(plot(Δλ,dotprods,seriestype=:scatter))
    readline()



    

    #print(size(evals))

    #eps = 10.0 .^ range(-9,-3)
    # eps = [1e-9]
    # mus = collect(-3:0.01:3.0)
    # ns = zeros(size(eps,1),size(mus,1))

    # for (l,e) in enumerate(eps)
    #     println(e)
    #     for (k, mu) in enumerate(mus)
    #         ns[l,k], mus[k] = count_evals(A,mu,e)
    #     end
    # end
    # display(plot(mus,[ns[l,:] for l=1:size(eps,1)],label=[l for l=1:size(eps,1)]))
    # readline() #prevents plot window from closing immediately (press enter in CLI to close plot window and terminate script)
end
