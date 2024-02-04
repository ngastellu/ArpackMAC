module CheckBandGaps

    include("./SpectralLanczos.jl")
    
    using .SpectralLanczos,  PyCall, Base.Filesystem, SparseArrays

    py"""
    import numpy as np
    vir_lbls = np.load('vir_lbls.npy')
    """

    lbls = PyArray(py"vir_lbls"o)
    rCC = 1.8

    ii = parse(Int, ARGS[1])
    n = lbls[ii]
  
    
    println("**** $n ****")
    eocc_file = "energies/occupied/eARPACK_bigMAC-$n.npy"
    eoccs = sort(PyArray(py"np.load($eocc_file)"o))

    evirt_file = "energies/virtual/eARPACK_bigMAC-$n.npy"
    evirts = sort(PyArray(py"np.load($evirt_file)"o))

    natoms_file = "natoms.npy"
    natoms = PyArray(py"np.load($natoms_file)"o)

    ii_file = "hamiltonians/inds/ii-$n.npy"
    iii = PyArray(py"np.load($ii_file)"o)

    jj_file = "hamiltonians/inds/jj-$n.npy"
    jj = PyArray(py"np.load($jj_file)"o)

    hvals_file = "hamiltonians/hvals/hvals-$n.npy"
    hvals = PyArray(py"np.load($hvals_file)"o)

    N = natoms[ii]
 
    println("Constructing hamiltonian...")
    H = sparse(iii,jj,hvals,N,N)
    H += H'
    println("Done!")


    if N % 2 == 0

        nhalf = Int(floor(N/2))

        k_occ = 0
        k_virt = 1

        eocc = eoccs[end-k_occ]
        evirt = evirts[k_virt]

        occ_Δn = zeros(Int, size(eoccs,1))
        virt_Δn = zeros(Int, size(evirts,1),)

        δn_occ, _ = count_evals(H, eocc)
        δn_virt, _ = count_evals(H, evirt)

        if δn_occ < nhalf - 1 && δn_vir > nhalf + 1
            println("Bad structure:")
            println("N/2 = $nhalf")
            println("δn_occ = $δn_occ")
            println("δn_virt = $δn_virt")

        else
            occ_Δn[k_occ+1] = δn_occ
            virt_Δn[k_virt] = δn_virt

            if δn_occ > nhalf - 1
                k_occ +=1
                eocc = eoccs[end-k_occ]
                δn_occ , _= count_evals(H,eocc)
                occ_Δn[k_occ+1] = δn_occ
            end

            if δn_virt < nhalf 
                k_virt +=1
                evirt = evirts[k_virt]
                δn_virt, _ = count_evals(H,evirt)
                virt_Δn[k_virt] = δn_virt
            end


            py"""
            nn = $n
            np.save(f'gap_check/occupied/odN-{nn}.npy', $(PyObject(occ_Δn)))
            np.save(f'gap_check/virtual/vdN-{nn}.npy', $(PyObject(virt_Δn)))
            """
        end

    else

        println("Odd number of atoms: N = $N")
    end

end