module Gershgorin

    using LinearAlgebra, SparseArrays

    export spectral_bounds

    function spectral_bounds(A; return_indiv_bounds=false)
        N = size(A,1)
        λmin = Inf
        λmax = -Inf

        if return_indiv_bounds
            indiv_bounds = zeros(2,N)
        end
        
        for i=1:N
            x0 = A[i,i]
            R = sum(abs.(A[i,:])) - x0
            lower = x0 - R
            upper = x0 + R
            if return_indiv_bounds
                indiv_bounds[:,i] = [lower; upper]
            end
            if lower < λmin
                λmin = lower
            end

            if upper > λmax
                λmax = upper
            end

        end
        if return_indiv_bounds
            return λmin, λmax, indiv_bounds
        else
            return λmin, λmax
        end
    end

    #function approx_dos(indiv_bounds)
    #end

end
