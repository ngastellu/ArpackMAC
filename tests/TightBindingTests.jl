module TightBindingTests

    include("../CoordsIO.jl")
    include("../TightBinding.jl")

    using Test, NPZ, .TightBinding, .CoordsIO

    istruc = 374
    pospath = "/Users/nico/Desktop/simulation_outputs/MAC_structures/kMC/slurm-6727121_fixed/sample-$(istruc).xsf"


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

    function loadpos() 
        pos, supercell = read_xsf(pospath)
        pos = pos[:,1:2]'
        for col in eachcol(pos)
            col .= col .% [40,40] # use `.=` to update cols in-place and not bind a new array to `col` (as would be the case with `=`, which turns this line into an assignment statement)
        end 
        return pos, supercell # transpose pos to have it column-ordered for more efficient looping
    end

    function load_nn_benchmarks()
        bench_data_dir = "nn_benchmarks/"
        dists = npzread(joinpath(bench_data_dir,"dists-$(istruc).npy"))
        nn_list = transpose(npzread(joinpath(bench_data_dir,"nns-$(istruc).npy")) .+ 1)
        dists_pbc = npzread(joinpath(bench_data_dir,"dists_pbc-$(istruc).npy"))
        nn_list_pbc = transpose(npzread(joinpath(bench_data_dir,"nns_pbc-$(istruc).npy")) .+ 1)

        return dists, nn_list, dists_pbc, nn_list_pbc
    end



    function nn_pairdists_vec_tst(N)
        pos, supercell = loadpos()
        dists, nn_list, dists_pbc, nn_list_pbc = load_nn_benchmarks()
        rCC = 1.8
        dists_new, ii, jj = nn_pairdists_vec(pos, rCC)
        nn_list_new = quick_hcat(ii, jj)

        println(nn_list_new[:,1:10])
        println(nn_list[:,1:10])

        println(dists[1:10])
        println(dists_new[1:10])
        println(size(dists))
        println(size(dists_new))
        # println(maximum(abs.(dists .- dists_new)))


        dists_pbc_new, ii, jj = nn_pairdists_vec(pos, rCC; cellsize=supercell)
        nn_list_pbc_new = quick_hcat(ii, jj)
        @test all(dists == dists_new)
        @test all(nn_list == nn_list_new)

        @test all(dists_pbc == dists_pbc_new)
        @test all(nn_list_pbc == nn_list_pbc_new)   
    end

    @testset "Nearest Neighbour list and distances test" begin
        nn_pairdists_vec_tst(istruc)
    end

end