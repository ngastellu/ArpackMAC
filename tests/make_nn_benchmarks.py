#!/usr/bin/env python

import numpy as np
import os
from scipy.spatial import KDTree
from time import perf_counter
from qcnico.coords_io import read_xsf, read_xyz

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

def get_supercell_ase(xyzpath):
    """Extracts supercell lattice dimensions from an XYZ file written by ASE.
        Assumes orthorhombic supercell; i.e. only extracts x-coord of 1st vector, y-coords of 2nd vector, z-coords of 3rd vector."""
    with open(xyzpath) as fo:
        fo.readline()
        line2 = fo.readline()
    supercell_coords = line2.split('=')[1].split() #get rid of trailing " characters
    orthorhombic_inds = [0,4,8]
    supercell = [float(supercell_coords[i].strip('"')) for i in orthorhombic_inds]
    return supercell


system = 'GNRs'

savedir = f"nn_benchmarks/{system}"
if not os.path.exists(savedir):
    os.makedirs(savedir)

if system == 'kMC_MAC':
    strucdir = "/Users/nico/Desktop/simulation_outputs/MAC_structures/kMC/slurm-6727121_fixed/"
    strucfiles = os.listdir(strucdir)

elif system == 'GNRs':
    strucdir = "/Users/nico/Desktop/simulation_outputs/GNRs/"
    lbls = ["armchair_11x50", "zigzag_11x100"]
    strucfiles = [f"gnr_{l}.xyz" for l in lbls]
else:
    raise ValueError(f"Invalid system type: {system}. Valid values are: ['kMC_MAC', 'GNRs'].")


savedir = f"nn_benchmarks/{system}"
if not os.path.exists(savedir):
    os.makedirs(savedir)

rCC = 1.8
k = 0

for sf in strucfiles:
    start = perf_counter()
    print(f'Working on file {sf}...', end = ' ')
    
    if sf.split('.')[-1] == 'xsf':
        lbl = sf.split('-')[-1].split('.')[0]
        pos, supercell = read_xsf(os.path.join(strucdir, sf),read_forces = False)
        supercell = np.array(supercell[:2])
        pos = pos[:,:2] % supercell
    else:
        lbl = "_".join(sf.split('.')[0].split("_")[1:])
        pospath = os.path.join(strucdir,sf) 
        pos = read_xyz(pospath)
        pos = pos[:,:2]
        supercell = get_supercell_ase(pospath)
        supercell = supercell[1:]

    np.save(os.path.join(savedir,f'pos-{lbl}.npy'), pos)

    read_end = perf_counter()
    read_time = read_end - start

    nn_inds, dists = nn_pairdists(pos,rCC)
    
    if k == 0:
        print("\nnn_inds = ", nn_inds)
        print("\ndists = ", dists)
        print(nn_inds.shape)
        print(dists.shape)

    np.save(os.path.join(savedir, f"dists-{lbl}.npy"),dists)
    np.save(os.path.join(savedir, f"nns-{lbl}.npy"),nn_inds)
    open_end = perf_counter()
    open_time = open_end - read_end

    print(pos.shape)
    print(supercell)

    nn_inds_pbc, dists_pbc = nn_pairdists(pos,rCC,cellsize=supercell)
    np.save(os.path.join(savedir, f"dists_pbc-{lbl}.npy"),dists_pbc)
    np.save(os.path.join(savedir, f"nns_pbc-{lbl}.npy"),nn_inds_pbc)
    pbc_end = perf_counter()
    pbc_time = pbc_end - open_end

    print('Done!')
    print(f'\tRead time = {read_time} seconds')
    print(f'\tNN dists time (no PBC) = {open_time} seconds')
    print(f'\tNN dists time (with PBC) = {pbc_time} seconds')

    k+=1