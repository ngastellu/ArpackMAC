module SpectralLanczos
  using LinearAlgebra, Arpack, SparseArrays, Printf, Statistics

  export spectral_shift_lanczos, count_evals, get_resids, sort_eigenpairs!, safe_eigs


    function safe_eigs(A::SparseMatrixCSC; nev=30, which=:LM, sigma=nothing,tol=0.0,maxiter=300,ritzvec=true,check=0,eps_shift=1e-8)
      # Wrapper for Arpack's `eigs` function which catches ZeroPivotError, if `sigma` is too close an eigenvalues
      try
        out = eigs(A;nev=nev,which=which,sigma=sigma,tol=tol,maxiter=maxiter,ritzvec=ritzvec,check=check)
      catch LoadError
        fail = true
        while fail
          if which == :LR
            println("[safe_eigs] ZPE encountered! old sigma = $sigma ---> new sigma =$(sigma+eps_shift)")
            sigma += eps_shift #shift up if we want to eigvals greater than `sigma`
          else 
            println("[safe_eigs] ZPE encountered! old sigma = $sigma ---> new sigma =$(sigma-eps_shift)")
            sigma -= eps_shift #shift down otherwise
          end
          out = eigs(A;nev=nev,which=which,sigma=sigma,tol=tol,maxiter=maxiter,ritzvec=ritzvec,check=check)
        end
      end
      return out
    end

    function count_evals(A::SparseMatrixCSC, shift::Number, eps_shift::Number=1e-9)
      # Counts the number of eigenvalues of A which are < mu
      l = ldlt(A,shift=-shift,check=false)
      while !issuccess(l)
        println("[count_evals] ZPE encountered! old shift = $shift ---> new shift =$(shift+eps_shift)")
        shift += eps_shift #eps_shift = a very small nb to add to the shift in case of a ZeroPivotError
        l = ldlt(A,shift=-shift,check=false) #this should fix the issue 
      end
      n_evals = sum(diag(l) .< 0)
      return n_evals, shift #return shift in case it has been modified due to ZeroPivotError
    end

    function orthogonalize(v0::Vector, V::AbstractArray)
      for j=1:size(V,2)
        v = V[:,j]
        if all(v .== zero(eltype(V)))
          println("Breaking at j = ", j)
          break
        end
        vproj = ( dot(v,v0)/(norm(v)*norm(v0)) ) * v
        v0 -= vproj
      end
      return v0/norm(v0)
    end


    function sort_eigenpairs!(evals,evecs,beta)
      inds = sortperm(evals)
      evals = evals[inds]
      evecs = evecs[:,inds]
      λM = evals[end]
      # Remove any eigenpairs > beta
      i = 0  
      while λM > beta 
        i+=1
        λM = evals[end-i]
      end
      if i > 0
        N = size(evals,1)
        println("[sort_eigenpairs!] Final i = ", N-i)
      end
      return evals[1:end-i], evecs[:,1:end-i]
    end

    function remove_duplicate_eigenpairs!(evals,evecs,r0,eps_ortho, conv_evals, conv_evecs)
      ninit = size(evals,1)
      println("Nb. Ritz pairs going in = ", ninit)
      dotprods = vec(abs.(r0' * evecs))
      println("Minimum dotprod = $(minimum(dotprods))")
      mask = dotprods .> eps_ortho #already converged eigenvectors are orthogonal to r0
      anti_mask = map(!,mask)
      if any(anti_mask)
        println("$(sum(anti_mask)) duplicates found. Listing them now: ")
        for (dp, lam, v) ∈ zip(dotprods[anti_mask], evals[anti_mask], eachcol(evecs[:,anti_mask]))
          diffs = abs.(lam .- conv_evals)
          mind = argmin(diffs)
          dots = vec(v' * conv_evecs)
          println(size(dots))
          println(size(diffs))
          Mind = argmax(dots)
          println("\ndp = $(dp) ; Min. diff to conv evals = $(diffs[mind]); argmin = $mind")#; $(diffs[Mind])")
          #println("Max dot prod w conv evecs = $(dots[Mind]); argmax = $Mind; $(dots[mind])")
          println("Corresponding dot prod w conv evecs = $(dots[mind])")
        end
        #print("Minimum failing dot product =  $(minimum(dotprods[anti_mask])) ; ")
      end
      if any(mask)
        print("Maximum passing dot product =  $(maximum(dotprods[mask])) ; ")
      end
      evecs = evecs[:,mask]  
      evals = evals[mask]
      println("Nb. Ritz pairs going out = $(size(evals,1)); ==> got rid of $(ninit - size(evals,1)) pairs.")
      #return evals, evecs
    end

    function dupcheck(evals,evecs,eps_ortho=1e-10,eps_diff=1e-9)
      n = size(evals,1)
      dups1 = Set{Tuple{Int,Int}}()
      dups2 = Set{Tuple{Int,Int}}()
      for (k, λ) ∈ enumerate(evals)
        diffs = abs.(λ .- evals)
        bbl = diffs .< eps_diff
        dotprods = vec(abs.(evecs[:,k]' * evecs))
        bbc = dotprods .> eps_ortho
        bbc[k], bbl[k] = false, false
        match = all(bbc .== bbl)
        if sum(!match) > 1
          mismatch_inds = findall(bbc .!= bbl)
          println("\n**** k = $k ****")
          println("!! Mismatched bool arrays !!")
          println("Mismatched inds = ", mismatch_inds)
          println("bbl[inds] = $(bbl[mismatch_inds])")
          println("diffs[inds] = $(diffs[mismatch_inds])")
          println("bbc[inds] = $(bbc[mismatch_inds])")
          println("dotprods[inds] = $(dotprods[mismatch_inds])")
        end
        d1 = findall(bbl)
        d1 = d1[d1 .> k] #look only ahead to avoid duplicate duplicates (lol)

        d2 = findall(bbc)
        d2 = d2[d2 .> k] #look only ahead to avoid duplicate duplicates (lol)
        for (d,D) ∈ [(d1,dups1), (d2,dups2)]
          for l ∈ d
            push!(D, (k,l))
          end
        end
      end
      println("\n **** TOTAL DUPLICATES FOUND = ($(length(dups1)), $(length(dups2))) ****")
      return dups1, dups2
    end

    function rde!(evals,evecs,r0,eps_ortho,conv_evals,conv_evecs)
      row1 = vec(abs.(evecs[1,:]))
      mask = row1 .> eps_ortho #already converged eigenvectors are orthogonal to r0
      anti_mask = map(!,mask)
      if any(anti_mask)
        println("$(sum(anti_mask)) duplicates found. Listing them now: ")
        for (dp, lam, v) in zip(row1[anti_mask], evals[anti_mask], eachcol(evecs[:,anti_mask]))
          diffs = abs.(lam .- conv_evals)
          mind = argmin(diffs)
          dots = vec(v' * conv_evecs)
          println(size(dots))
          println(size(diffs))
          Mind = argmax(dots)
          println("\ns1i = $(dp) ; lambda $(lam) ; Min. diff to conv evals = $(diffs[mind]); argmin = $mind")#; $(diffs[Mind])")
          #println("Max dot prod w conv evecs = $(dots[Mind]); argmax = $Mind; $(dots[mind])")
          println("Corresponding dot prod w conv evecs = $(dots[mind])")
        end
        print("Minimum failing s1i =  $(minimum(row1[anti_mask])) ; ")
      end
      if any(mask)
        print("Minimum s1i =  $(minimum(row1[mask])) ; ")
      end
      evecs = evecs[:,mask]  
      evals = evals[mask]
    end

    function rde2!(evals, evecs, eps_ortho, conv_evals, conv_evecs)
      dotprods = evecs' * conv_evecs
      mask = dotprods .< eps_ortho #already converged eigenvectors are orthogonal to r0
      anti_mask = map(!,mask)
      if any(anti_mask)
        println("Duplicates found. Listing them now: ")
        for m=1:size(evals,1)
          bbb = vec(anti_mask[m,:])
          if any(bbb)
            lam = evals[m]
            diffs = abs.(lam .- conv_evals)
            diffinds = sortperm(diffs)
            print("\nlambda = $lam: ")
            println(" Nb. of possible duplicates = $(sum(bbb))")
            println("Possible dups: $(vec(conv_evals[bbb]))")
            println("Dot prods = $(dotprods[m,bbb])")
          end
        end 
      end
      mask = vec(reduce(&,mask,dims=2))
      evecs = evecs[:,mask]  
      evals = evals[mask]
    end

    function rdef(evecs,evals,conv_evecs,eps_ortho)
      for V ∈ [evecs, conv_evecs]
        dotprods = abs.(evecs' * V)
        if V === evecs
          dotprods[diagind(dotprods)] .= 0 # ignore dot product w/ self
        else
          eps_ortho = 0.01 #orthogonality condition needs to be STRONGLY relaxed when comparing Ritz vecs from different iterations
        end
        good_mask = vec(all(dotprods .< eps_ortho, dims=2)) #duplicate eigenvector is one whose dot product is nonzero with other eigenvectors
        evecs = evecs[:,good_mask]
        evals = evals[good_mask]
      end
      return evals, evecs
    end

    function rdef_debug(evals,evecs,conv_evals,conv_evecs,eps_ortho)
      ninit = size(evals,1)
      println("Nb. Ritz pairs going in = ", ninit)
      for V ∈ [evecs, conv_evecs]
        dotprods = abs.(evecs' * V)
        if V === evecs
          println("~~~ Doin evecs ~~~")
          dotprods[diagind(dotprods)] .= 0 # ignore dot product w/ self
          Λ = evals
        else
          println("~~~ Doin conv_evecs now ~~~")
          Λ = conv_evals
          eps_ortho = 0.01 #orthogonality condition needs to be STRONGLY relaxed when comparing Ritz vecs from different iterations
        end
        dpbools = dotprods .< eps_ortho
        good_mask = vec(all(dpbools, dims=2)) #duplicate eigenvector is one whose dot product is nonzero with other eigenvectors
        if !all(good_mask)
          bad_mask = map(!,good_mask)
          for j ∈ findall(bad_mask)
            println("\n\n*** j = $j ***")
            badinds = vec(findall(map(!,dpbools[j,:])))
            println("Bad dotprods = ", dotprods[j,badinds])
            println("ediff = ", abs.(evals[j] .- Λ[badinds]))
          end
        end
        evecs = evecs[:,good_mask]
        evals = evals[good_mask]
      end
      println("Nb. Ritz pairs going out = $(size(evals,1)); ==> got rid of $(ninit - size(evals,1)) pairs.")
      return evals, evecs
    end

    function updateC!(conv_evals,conv_evecs,evals,evecs,nC)
      # This function updates the set C of converged eigenpairs with newly found eigenvalues/vectors
      # Pragmatically, C is really two sets: one containing the eigenvalues, the other containing the eigenvectors
      nconv = size(evals,1)
      conv_evecs[:,1+nC:nC+nconv] = evecs
      conv_evals[1+nC:nC+nconv] = evals
      nC += nconv
    end

    function updateC_debug(conv_evals,conv_evecs,evals,evecs,nC)
      # This function updates the set C of converged eigenpairs with newly found eigenvalues/vectors
      # Pragmatically, C is really two sets: one containing the eigenvalues, the other containing the eigenvectors
      nconv = size(evals,1)
      nC += nconv
      if nC > size(conv_evals,1)
        N = size(evecs,1)
        println("!!!! OVERFLOW !!!! Reallocating... converged arrays.")
        tmp_evals = fill(fill_value,nC)
        tmp_evecs = zeros(N,nC)

        tmp_evals[1:nC-nconv] = conv_evals[:]
        tmp_evals[nC-nconv:nC] = evals[:]

        tmp_evecs[:,1:nC-nconv] = conv_evecs[:,:]
        tmp_evecs[:,nC-nconv:nC] = evecs[:,:]

        conv_evals = tmp_evals
        conv_evecs = tmp_evecs
      end
      conv_evecs[:,1+nC:nC+nconv] = evecs
      conv_evals[1+nC:nC+nconv] = evals

      return conv_evals, conv_evecs, nC
    end

    function get_resids(evals,evecs,A)
      nconv = size(evals,1)
      u1 = A * evecs
      u2 = evecs * diagm(evals)
      resid = zeros(nconv)
      for j=1:nconv
        resid[j] = norm(u1[:,j] - u2[:,j])
      end
      return resid
    end
  
    function spectral_shift_lanczos(A::SparseMatrixCSC, alpha::Number, beta::Number, 
      nlanczos::Integer, eps_lanczos=1e-9, eps_ortho::Number=1e-14,eps_diff=1e-9)
        # Get target number of eigenvalues
        println("Getting ntarget...")
        nA,alpha = count_evals(A,alpha)
        nB,beta = count_evals(A,beta)
        ntarget = nB - nA
        println("Done! ntarget = ", ntarget)
        print('\n')

        #Initialise variables
        N = size(A,1)
        r0 = randn(N)
        r0 = r0/norm(r0)
        fill_value = 1000*beta
        conv_evals = fill(float(fill_value), ntarget) #initialise evals with 1000*β ∉[α,β] to avoid collisions down the line 
        conv_evecs = zeros(N,ntarget)
        nC = 0
        r = Int(nlanczos/2.5) - 2 #number of deisred eigenpairs per step
        println("r = ", r)

        # Get first batch of eigenvalues
        print("Getting 1st batch of eigenpairs... ")
        evals, evecs, nconv = eigs(A,nev=Int(r/2),sigma=alpha,
          maxiter=nlanczos,which=:LR,v0=r0,tol=eps_lanczos,check=2)[1:3]
        print("Done! Sorting... ")
        sort_eigenpairs!(evals,evecs,beta)
        print("Done! Updating C... ")
        nC = updateC!(conv_evals,conv_evecs,evals,evecs,nC)
        println("Done!")



        # *********** Begin iteration accross [alpha;beta] ***********
        println("\nCommencing loop.")
        mu_prev = alpha
        mu = mu_prev + 2*(evals[end] - mu_prev)
        nleft_prev = nA
        o=1
        while nC < ntarget #&& o < 3
          println("\n\n----------------------- o = $o -----------------------")
          o+=1
          println("\nmu = $mu ; mu_prev = $mu_prev ; nC = $nc\n")
          nleft, mu = count_evals(A, mu)
          ntarget_mu = nleft - nleft_prev #number of eigenvalues between mu_prev and mu
          println("ntarget_mu = ", ntarget_mu)
          #println("ntarget = ")

          # Redefine shift in case too many eigenvalues lie in [mu_prev;mu]
          while ntarget_mu > 10 * r
            mu = mu - 0.5*(mu - mu_prev) #dichotomise [mu_prev, mu]
            nleft, mu = count_evals(A, mu)
            ntarget_mu = nleft - nleft_prev
          end
          
          # Get new random starting vector, orthogonalize with already-converged eigenvectors
          #print("Getting new r0... ")
          r0 = orthogonalize(randn(N),conv_evecs)
          #println("Done!")

          # Run Lanczos to get new eigenpairs
          evals, evecs, nconv = eigs(A,nev=r,maxiter=nlanczos,sigma=mu,
            which=:LM,v0=r0,tol=eps_lanczos,check=2)[1:3]
          sort_eigenpairs!(evals,evecs,beta)
          print("Removing duplicates... ")
          #remove_duplicate_eigenpairs!(evals,evecs,r0,eps_ortho,conv_evals,conv_evecs)
          evals, evecs = rdef_debug(evals,evecs,conv_evals,conv_evecs,eps_ortho)
          #rdef!(evecs,conv_evecs,eps_ortho)
          println("Done!")
          println("Size of evals now = $(size(evals,1))")
          #println("Checking duplicates...")
          #dupcheck(evals,evecs,eps_diff)
          #println("Beta - max(eval) = ", beta - evals[end])
          #println("Any evals > beta: ", any(evals .> beta))
          nC = updateC!(conv_evals,conv_evecs,evals,evecs,nC)
          #conv_evals, conv_evecs, nC = updateC_debug(conv_evals,conv_evecs,evals,evecs,nC)

          right_check = conv_evals .> mu #check if we have any prev eigenvalues > mu
          left_check = conv_evals .< mu #idem for eigenvalues < mu
          println("Nb. of evals < mu = ", sum(evals .< mu))
          println("Nb. of evals > mu = ", sum(evals .> mu))
          println("Nb. of evals == mu = ", sum(evals .== mu))

          if any(right_check) && any(left_check) #ideal case
            # Get eigenpairs between mu_prev and mu
            mask_mu = (conv_evals .< mu) .* (conv_evals .> mu_prev)
            evals_mu = @view conv_evals[mask_mu]
            evecs_mu = @view conv_evecs[:,mask_mu]
            nconv_mu = size(evals_mu,1)
            println("nconv_mu = ", nconv_mu)
            delta = ntarget_mu - nconv_mu #nb of eigenpairs left to obtain in [mu_prev; mu]
            println("delta = ", delta)
            if delta < 0
              println("****** ALARM: nconv_mu > ntarget_mu ******")
              #return conv_evals, mu, mu_prev
            end
            while delta > 0
              println("Missing evals in [mu_prev;mu]. Re-running Lanczos w/ new |r0>.")
              print("Getting new |r0>... ")
              r0 = orthogonalize(randn(N),evecs_mu)
              print("Done! ")
              print("Running Lanczos... ")
              evals2, evecs2 = eigs(A,nev=r,maxiter=nlanczos,sigma=mu,
                which=:SR,v0=r0,tol=eps_lanczos,check=2)[1:3]
              print("Done! ")
              print("Removing duplicates... ")
              #remove_duplicate_eigenpairs!(evals2,evecs2,r0,eps_ortho,conv_evals,conv_evecs)
              #rdef!(evecs2,conv_evecs,eps_ortho)
              evals2, evecs2 = rdef_debug(evals2,evecs2,conv_evals,conv_evecs,eps_ortho)
              #ndups = dupcheck(evals2,evecs2)
              print("Done! ")
              nconv = size(evals2,1)
              println("nconv = ", nconv)
              #cat!(evecs, evecs2)
              #append!(evals,evals2)
              if nconv > 0
                sort_eigenpairs!(evals2,evecs2,beta)
                nC = updateC!(conv_evals,conv_evecs,evals2,evecs2,nC)
                #conv_evals, conv_evecs, nC = updateC_debug(conv_evals,conv_evecs,evals2,evecs2,nC)
              end
              delta -= nconv
            end
            mu_prev = mu
            nleft_prev = nleft
            mu = mu + 2 * (maximum(evals) - mu)
            #updateC!(conv_evals,conv_evecs,evals,evecs,nC)

          elseif any(right_check) && ntarget_mu > 0 #no converged eigenpairs < mu
            println("No eigenvalues < mu; and ntarget_mu > 0 ...")
            mu2 = mu_prev + 1/(mu - minimum(evals)) #define intermediate mu2 to assert mu2 < mu
            println("mu2 < mu: ", mu2 < mu)
            mu = mu2
            nleft = nleft_prev # reset nleft to previous value
            # Still add converged eigenpairs to C 
            #nC = updateC!(conv_evals,conv_evecs,evals,evecs,nC)

          elseif any(left_check) && evals[end] < beta && mu < beta # no eigenpairs converged > mu
            println("No eigenvalues > mu...")
            # Update converged eigenpairs
            #nC = updateC!(conv_evals,conv_evecs,evals,evecs,nC)
            # Reset tolerance to get under-converged eigenvalues, get new shift, and retry
            old_tol = eps_lanczos
            while !any(right_check)
              new_tol = old_tol * 100
              evals2 = eigs(A,nev=10*r,maxiter=nlanczos*10,sigma=mu,
                which=:LR,v0=r0,tol=new_tol,ritzvec=false,check=2)[1]
              right_check = evals2 .> mu 
            end
            mu_prev = mu
            nleft_prev = nleft
            mu = maximum(evals2)

          else
            println("~~~~~~~~~~ TWILIGHT ZONE ~~~~~~~~~~ (You shouldn't see this -- not good.)")
          end

        end
        println("\n\n----------------------- Finished main loop -----------------------")
        final_inds = conv_evals .!= fill_value
        conv_evals = conv_evals[final_inds]
        conv_evecs = conv_evecs[:,final_inds]
        #sort_eigenpairs!(conv_evals,conv_evecs,beta)
        println("Final dupcheck...")
        ndups = dupcheck(conv_evals,conv_evecs)
        return conv_evals, conv_evecs, alpha, beta, mu_prev #return alpha, beta in case of ZeroPivotError when running ldlt

      end
end