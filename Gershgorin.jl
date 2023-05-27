module Gershgorin

    using LinearAlgebra, SparseArrays

    export spectral_bounds

    function spectral_bounds(A)
        N = size(A,1)
        λmin = -Inf
        λmax = Inf
        indiv_bounds = zeros(2,N)
        for i=1:N
            x0 = A[i,i]
            R = sum(A[i,:]) - x0
            lower = x0 - R
            upper = x0 + R
            indiv_bounds[:,i] = [lower; upper]
            if lower > λmin
                λmin = lower
            end

            if upper < λmax
                λmax = upper
            end

        end
        return λmin, λmax, return indiv_bounds
    end

    #function approx_dos(indiv_bounds)
    #end

end
