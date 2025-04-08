#!/usr/bin/env python

import numpy as np
import os
from scipy.spatial import KDTree
from time import perf_counter
from qcnico.coords_io import read_xsf

def nn_pairdists(pos,rCC,cellsize=None):
    tree = KDTree(pos,boxsize=cellsize)
    nn_dict = tree.sparse_distance_matrix(tree, rCC)
    dists = np.array(list(nn_dict.values()))
    nn_inds = np.array(list(nn_dict.keys()))
    
    # keep only upper triangular entries (i.e. above main diagonal of distance matrix)
    diffs = np.diff(nn_inds, axis=1)
    upper_tri = (diffs > 0).flatten()
    nn_inds_unsorted = nn_inds[upper_tri]
    dists_unsorted = dists[upper_tri]

    # Now sort index pairs (i,j) stored in `nn_inds_unsorted` by i, then by j.
    # `lexsort`` sorts by LAST key first, so we use `roll` to sort by rows first
    sorted_inds = np.lexsort(np.roll(nn_inds_unsorted,1,axis=1).T) 

    return nn_inds_unsorted[sorted_inds], dists_unsorted[sorted_inds]


strucdir = "/Users/nico/Desktop/simulation_outputs/MAC_structures/kMC/slurm-6727121_fixed/"

savedir = "nn_benchmarks"
if not os.path.exists(savedir):
    os.makedirs(savedir)

rCC = 1.8
k = 0

for xsf in os.listdir(strucdir):
    
    start = perf_counter()
    print(f'Working on file {xsf}...', end = ' ')
    istruc = xsf.split('-')[-1].split('.')[0]
    pos, supercell = read_xsf(os.path.join(strucdir, xsf),read_forces = False)
    supercell = np.array(supercell[:2])
    pos = pos[:,:2] % supercell
    read_end = perf_counter()
    read_time = read_end - start

    nn_inds, dists = nn_pairdists(pos,rCC)
    
    if k == 0:
        print("\nnn_inds = ", nn_inds)
        print("\ndists = ", dists)
        print(nn_inds.shape)
        print(dists.shape)

    np.save(os.path.join(savedir, f"dists-{istruc}.npy"),dists)
    np.save(os.path.join(savedir, f"nns-{istruc}.npy"),nn_inds)
    open_end = perf_counter()
    open_time = open_end - read_end

    print(pos.shape)
    print(supercell)

    nn_inds_pbc, dists_pbc = nn_pairdists(pos,rCC,cellsize=supercell)
    np.save(os.path.join(savedir, f"dists_pbc-{istruc}.npy"),dists_pbc)
    np.save(os.path.join(savedir, f"nns_pbc-{istruc}.npy"),nn_inds_pbc)
    pbc_end = perf_counter()
    pbc_time = pbc_end - open_end

    print('Done!')
    print(f'\tRead time = {read_time} seconds')
    print(f'\tNN dists time (no PBC) = {open_time} seconds')
    print(f'\tNN dists time (with PBC) = {pbc_time} seconds')

    