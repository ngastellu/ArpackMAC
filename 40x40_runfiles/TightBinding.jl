module TightBinding

using SparseArrays, LinearAlgebra

export constructHtb_dense, lindbergHtb_sparse, nn_pairdists, nn_pairdists_vec, nn_pairdists_full_vec

function simpleHtb_dense(n::Int64,t::Number)
    H = zeros((n,n))
    for j=1:n-1
        H[j,j+1] = t
        H[j+1,j] = t
    end
    H[end,1] = t
    H[1,end] = t
    return H
end

constructHtb_dense = simpleHtb_dense

function ij_inds(k_inds,N)
    # This function returns the indices (i,j) labelling a pair of atoms
    # such that pairdists[k] = norm(pos[i]-pos[j]), where pairdists is a 
    # 1D array constructed in the `nn_pairdists` function.
    zero_j_inds = [Int(k*(k-1)/2) for k=2:N]
    i_inds = [sum(k .>= zero_j_inds) + 1 for k ∈ k_inds]
    j_inds = @. k_inds - zero_j_inds[i_inds-1] + 1
    return i_inds, j_inds
end

function ij_inds_vec(k_inds,N)
    # This function returns the indices (i,j) labelling a pair of atoms
    # such that pairdists[k] = norm(pos[i]-pos[j]), where pairdists is a 
    # 1D array constructed in the `nn_pairdists` function.
    N_inds = [n*N - Int((n*(n+1)/2)) for n=1:N-1]
    ii = zeros(Int,size(k_inds,1))
    jj = zeros(Int,size(k_inds,1))
    n = 1
    for (m, k) in enumerate(k_inds)
        if k > N_inds[n]
            n+=1
        end
        ii[m] = n 
        jj[m] = k + N - N_inds[n]
    end
    # ii = [sum(k .>= N_inds) + 1 for k ∈ k_inds]
    # jj = @. Int(k_inds + (1 - ii)*N + ii*(ii-1)/2) + 1
    return ii, jj
end


function nn_pairdists(pos, rNN)
    N = size(pos,1)
    pairdists = zeros(Int(N*(N-1)/2)) 
    k = 1
    for i=2:N 
        for j=1:i-1
            pairdists[k] = norm(pos[i,:]-pos[j,:])
            k+=1
        end
    end
    bonded = pairdists .<= rNN
    nn_k_inds = findall(bonded)
    ii, jj = ij_inds(nn_k_inds, N)
    return pairdists[bonded], ii, jj
end

function nn_pairdists_vec_ht(pos, rNN, return_hashtable=true)
    N = size(pos,1)
    pairdists = zeros(Int(N*(N-1)/2))
    if return_hashtable
        hashtable = zeros(Int,Int(N*(N-1)/2),2)
    end
    k = 1
    for i=1:N-1 
        pairdists[k:k+N-i-1] = norm.(eachrow(pos[i+1:N,:] .- pos[i,:]'))
        if return_hashtable
            hashtable[k:k+N-i-1,:] = reduce(hcat,[[i,j] for j=i+1:N])'
            
        end
        k+=N-i
    end
    bonded = pairdists .<= rNN
    if return_hashtable
        #hashtable = reduce(hcat, hashtable)' #turns hashtable into a matrix (as opposed to a vector of vectors)
        return pairdists[bonded], hashtable[bonded,:]
    else
        nn_k_inds = findall(bonded)
        ii, jj = ij_inds_vec(nn_k_inds, N)
        return pairdists[bonded], ii, jj
    end
end

function nn_pairdists_vec(pos, rNN)
    N = size(pos,1)
    pairdists = zeros(Int(N*(N-1)/2))
    ii = zeros(Int,Int(N*(N-1)/2))
    jj = zeros(Int,Int(N*(N-1)/2))
    k = 1
    for i=1:N-1 
        pairdists[k:(k+N-i-1)] = norm.(eachrow(pos[i+1:N,:] .- pos[i,:]'))
        ii[k:(k+N-i-1)] .= i 
        jj[k:(k+N-i-1)] = i+1:N
        k+=N-i
    end
    bonded = pairdists .<= rNN
    return @views pairdists[bonded], ii[bonded], jj[bonded]
end

function nn_pairdists_full_vec(pos,rNN)
    N = size(pos,1)
    pairdists = norm.(pos' .- pos, 2)
    for i=1:N
        pairdists[i,i] = rNN * 10 #ignore diagonal elements when determining which pairs are NN 
    end
    bonded = findall(pairdists .<= rNN)
    nn_pairdists = pairdists[bonded]
    ii = [ci[1] for ci in bonded]
    jj = [ci[2] for ci in bonded]
    return nn_pairdists, ii, jj
end

function lindbergHtb_sparse(pos,rNN)
    β0 = -2.438 #eV
    kb = 0.405 #angstrom^-1
    R0 = 1.397 #angstrom
    μb = 2.035 #angstrom^-1

    N = size(pos,1)
    println("Entering pairdists now...")
    dists, ii, jj = nn_pairdists_vec(pos,rNN)
    hvals = @. β0 * exp(-μb*(dists-R0)) * (1+kb*(dists-R0))
    Htb = sparse(ii,jj,hvals,N,N)
    Htb += Htb' #symmetrise Htb
end

end
