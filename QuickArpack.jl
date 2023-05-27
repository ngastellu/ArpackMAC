module QuickArpack
include("./TightBinding.jl")
include("./SpectralLanczos.jl")
using .TightBinding, .SpectralLanczos
using LinearAlgebra, SparseArrays, Arpack, Plots

N = 10000
t = -1
H = sparse(constructHtb_dense(N,t))

nmax = 500
nstep=50
nrange = N .- collect(0:nstep:nmax)
nsamples = size(nrange,1)

rdat = fill(8000.0, (nsamples,3)) 
nconvs = fill(8000,nsamples)
nevals = fill(8000,nsamples)
niters = fill(8000,nsamples)

for (k,n) ∈ enumerate(nrange)
    if n == 0
        n = 1
    end
    print("$n :  ")
    e,ϕ,nconv,niter, _, _ = eigs(H,nev=n,sigma=1e-9,which=:LR,tol=1e-9,maxiter=1000000,check=1)
    sort_eigenpairs!(e,ϕ,1000)
    α = e[1] - 1e-9 #shift boundaries to avoid ZPEs
    β = e[end] + 1e-9 
    nα, α = count_evals(H,α)
    nβ, β = count_evals(H,β)
    nevals[k] =nβ - nα # 'true' number of eigenvalues between the largest and smallest eigenvalues that have converged
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
p3 = plot(nrange,nevals .- nconvs, label="# missing eigenpairs")
p = plot(p1,p2,p3,layout=(3,1))
savefig(p,"arpack_test.pdf")
display(p)
readline()

n = 1:N
realE = @. -2*cos(2*π*n/N)

display(plot(n,realE))
readline()

end