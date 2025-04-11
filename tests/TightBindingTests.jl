module TightBindingTests

    include("../CoordsIO.jl")
    include("../TightBinding.jl")

    using Test, NPZ, .TightBinding, .CoordsIO

    # Define useful globals
    rCC = 1.8

    # MAC globals
    istruc = 374
    pospath_MAC = "/Users/nico/Desktop/simulation_outputs/MAC_structures/kMC/slurm-6727121_fixed/sample-$(istruc).xsf"
    bench_data_dir_MAC = "nn_benchmarks/kMC_MAC/"

    # GNR globals        
    gnr_label = "zigzag_11x100"
    gnr_type = split(gnr_label, '_')[1]
    gnr_pospath = "/Users/nico/Desktop/simulation_outputs/graphene_TB/$(gnr_type)/gnr_$(gnr_label).xyz"
    bench_data_dir_gnr = "nn_benchmarks/GNRs/"
    
    function quick_hcat(ii, jj)
        N = size(ii, 1)
        @assert size(jj, 1) == N
        nn = Array{Int64}(undef,2,N)

        for k=1:N
            nn[1,k] = ii[k]
            nn[2,k] = jj[k]
        end
        return nn
    end

    function loadpos_MAC(pospath) 
        pos, supercell = read_xsf(pospath)
        pos = pos[:,1:2]'
        supercell = supercell[1:2]
        for col in eachcol(pos)
            col .= mod.(col, supercell) # use `.=` to update cols in-place and not bind a new array to `col` (as would be the case with `=`, which turns this line into an assignment statement)
        end 
        return pos, supercell # transpose pos to have it column-ordered for more efficient looping
    end

    function loadpos_gnr(pospath)
        pos, supercell, _ = read_xyz_supercell(pospath)
        pos = pos[1:2,:]
        supercell = supercell[2:3]
        return pos, supercell
    end

    function load_nn_benchmarks(bench_data_dir, label)
        dists = npzread(joinpath(bench_data_dir,"dists-$(label).npy"))
        nn_list = transpose(npzread(joinpath(bench_data_dir,"nns-$(label).npy")) .+ 1)
        dists_pbc = npzread(joinpath(bench_data_dir,"dists_pbc-$(label).npy"))
        nn_list_pbc = transpose(npzread(joinpath(bench_data_dir,"nns_pbc-$(label).npy")) .+ 1)

        return dists, nn_list, dists_pbc, nn_list_pbc
    end



    function nn_pairdists_vec_tst(pos, supercell, bench_data_dir, label)
        dists, nn_list, dists_pbc, nn_list_pbc = load_nn_benchmarks(bench_data_dir, label)
        
        dists_new, ii, jj = nn_pairdists_vec(pos, rCC)
        nn_list_new = quick_hcat(ii, jj)

        dists_pbc_new, ii, jj = nn_pairdists_vec(pos, rCC, supercell)
        nn_list_pbc_new = quick_hcat(ii, jj)

        println(nn_list_pbc_new[:,1:10])
        println(nn_list_pbc[:,1:10])

        println(dists_pbc[1:10])
        println(dists_pbc_new[1:10])
        println(size(dists_pbc))
        println(size(dists_pbc_new))

        @test all(dists .== dists_new)
        @test all(nn_list .== nn_list_new)

        @test all(dists_pbc .== dists_pbc_new)
        @test all(nn_list_pbc .== nn_list_pbc_new)   
    end

    @testset "Positions test" begin
        pos, supercell = loadpos_MAC(pospath_MAC)
        pos_npy = (npzread(joinpath(bench_data_dir_MAC, "pos-$(istruc).npy")))'
        @test all(pos == pos_npy)
    end

    @testset "[MAC] NN list and distances test" begin
        pos, supercell = loadpos_MAC(pospath_MAC)
        nn_pairdists_vec_tst(pos, supercell, bench_data_dir_MAC,istruc)
    end

    @testset "[GNR] NN list and distances test" begin
        pos, supercell = loadpos_gnr(gnr_pospath)
        nn_pairdists_vec_tst(pos, supercell, bench_data_dir_gnr, gnr_label)
    end

end