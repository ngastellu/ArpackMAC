#!/usr/bin/env python

import numpy as np
from qcnico.coords_io import read_xsf
from scipy.spatial import cKDTree

def gen_perm(N,n,m):
    n = np.array(n)
    m = np.array(m)
    P = np.eye(N,dtype=int)
    P[n,n] = 0
    P[m,m] = 0
    P[n,m] = 1
    P[m,n] = 1
    return P

posfile = "data/bigMAC_10x10-64_relaxed.xsf"
pos, _ = read_xsf(posfile,read_forces=False)
N = pos.shape[0]
rCC = 1.8
tree = cKDTree(pos)
nn_pairs = np.sort(np.vstack(list(tree.query_pairs(rCC))), axis=1)
print(nn_pairs.shape)
ii = np.argsort(nn_pairs[:,0])
nn_pairs = nn_pairs[ii,:]
np.save('nn_inds.npy',nn_pairs)
diffs = np.diff(nn_pairs,axis=0)
unsortedbools = ((diffs[:,0] == 0) * (diffs[:,1] < 0)).nonzero()[0]
P = gen_perm(nn_pairs.shape[0],unsortedbools, unsortedbools+1)
print(np.all(P.sum(1) == 1))
nn_pairs = P @ nn_pairs
print(nn_pairs.dtype)
print(np.all(np.diff(nn_pairs,axis=0)>=0))
dists = np.linalg.norm(pos[None,:] - pos[:,None], axis=2)
np.save('nn_dists.npy',dists[nn_pairs[:,0], nn_pairs[:,1]])