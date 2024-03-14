#!/usr/bin/env python

import numpy as np 
from qcnico.coords_io import read_xsf, read_xyz
from qcnico.remove_dangling_carbons import remove_dangling_carbons
from scipy.spatial import cKDTree


strucdir = '/Users/nico/Desktop/simulation_outputs/MAC_structures/Ata_structures/t1/'
strucdir_relaxed = strucdir + 'relaxed/'


nstrucs = 38
rCC = 1.8
eps_dist = 1e-10

np.random.seed(0)
istrucs = np.random.choice(nstrucs,size=10,replace=False)

for n in istrucs:
    print(n)
    pos = read_xyz(strucdir + f't1n{n}.xyz')
    pos_relaxed, _ = read_xsf(strucdir_relaxed + f't1n{n}_relaxed.xsf')

    pos = remove_dangling_carbons(pos,rCC)
    pos_relaxed = remove_dangling_carbons(pos_relaxed,rCC)

    # assert pos.shape == pos_relaxed.shape, 'Shapes don\'t match!'
    tr = cKDTree(pos)
    tr_relaxed = cKDTree(pos_relaxed)
    
    v_close_atoms = tr.query_ball_tree(tr_relaxed, eps_dist)
    nequivalent_atoms = np.array([len(l) for l in v_close_atoms])
    for k, iatoms in enumerate(v_close_atoms):
        if len(iatoms) > 0:
            print(f'{k} --> {iatoms}')
    print('Each atom has an equivalent relaxed image: ', np.all(nequivalent_atoms == 1))
    


    """ This brute force approach is too computationally expensive... woops. """
    # dists = np.linalg.norm(pos[:,None,:] - pos_relaxed,axis=2)
    # izero_dists = (dists==0).nonzero()
    # print(len(izero_dists))
    # print(pos.shape[0])
