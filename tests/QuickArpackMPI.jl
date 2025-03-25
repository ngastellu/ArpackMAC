module QuickArpackMPI
include("./TightBinding.jl")
include("./SpectralLanczos.jl")
using .TightBinding, .SpectralLanczos
using LinearAlgebra, SparseArrays, Arpack, PyCall,  MPI

MPI.Init()
comm = MPI.COMM_WORLD

rank = MPI.Comm_rank(comm) + 1
nprocs = MPI.Comm_size(comm)

# NOT DONE HERE

nmax = 2000
nstep=50
nrange = collect(0:nstep:nmax)
nsamples = size(nrange,1)

rdat = fill(8000.0, (nsamples,3)) 
nconvs = fill(8000,nsamples)
niters = fill(8000,nsamples)

N = 10000
t = -1
H = sparse(constructHtb_dense(N,t))


for (k,n) ∈ enumerate(nrange)
    if n == 0
        n = 1
    end
    print("$n :  ")
    e,ϕ,nconv,niter, _, _ = eigs(H,nev=n,sigma=1e-9,which=:SR,tol=1e-9,maxiter=1000000,check=1)
    nconvs[k] = nconv
    niters[k] = niter
    resids = get_resids(e,ϕ,H)
    inds = sortperm(resids)
    rdat[k,1] = resids[inds[1]] #minimum residual
    rdat[k,2] = sum(resids)/nconv #average residual
    rdat[k,3] = resids[inds[end]] #max residual
    worst_e = e[inds[1]]
    best_e = e[inds[end]]
    println("best_e = $(best_e) ; worst_e = $(worst_e)")
end

p1 = plot(nrange,[rdat[:,k] for k=1:3],label=["Lowest" "Mean" "Highest"])
p2 = plot(nrange,niters, label="# iterations")
p3 = plot(nrange,nconvs, label="# converged")
p = plot(p1,p2,p3,layout=(3,1))
savefig(p,"arpack_test.pdf")
display(p)
readline()

end