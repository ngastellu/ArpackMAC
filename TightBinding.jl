module TightBinding

using SparseArrays, NearestNeighbors

export constructHtb_dense, lindbergHtb_sparse

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
    zero_j_inds = [Int(k*(k-1)/2) for k=1:N]
    i_inds = [sum(k .>= zero_j_inds) for k ∈ k_inds]
    j_inds = k .- zero_j_inds[i_inds .- 1]
    return i_inds, j_inds
end


function nn_pairdists(pos, rNN)
    N = size(pos,1)
    pairdists = zeros(Int(N*(N-1)/2))
    k = 1
    for i=1:N 
        for j=1:i-1
            pairdists[k] = norm(pos[i]-pos[j])
        end
    end
    bonded = pairdists .<= rNN
    nn_k_inds = findall(bonded)
    ii, jj = ij_inds(nn_k_inds)
    return pairdists[bonded], ii, jj
end

function lindbergHtb_sparse(pos,rNN,gamma)
    b0 = -2.438 #eV
    kb = 0.405 #angstrom^-1
    R0 = 1.397 #angstrom
    mub = 2.035 #angstrom^-1

    N = size(pos,1)
    dists, ii, jj = nn_pairdists(pos,rNN)
    @. hvals = b0 * exp(-mub*(dists-R0)) * (1+kb*(dists-R0))
    Htb = sparse(ii,jj,hvals)
end

end
