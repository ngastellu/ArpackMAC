#!/usr/bin/env python


from glob import glob
import numpy as np
from scipy.spatial import KDTree
from qcnico.coords_io import read_xyz, read_xsf
from qcnico.remove_dangling_carbons import remove_dangling_carbons



def get_nn_dists(posfile):
    rCC = 1.8
    if posfile.split('.')[-1] == 'xyz':
        pos = remove_dangling_carbons(read_xyz(posfile),rCC)
    elif posfile.split('.')[-1] == 'xsf':
        pos = remove_dangling_carbons(read_xsf(posfile)[0],rCC)
    tree = KDTree(pos)
    M = tree.sparse_distance_matrix(tree,rCC)
    nn_dists = M.nnz_elems # <---- fix this
    return nn_dists



strucdir = '/Users/nico/Desktop/simulation_outputs/MAC_structures/'

t1_unrel_strucs = glob(strucdir + 'Ata_structures/*xyz')
t1_rel_strucs = glob(strucdir + 'Ata_structures/relaxed/*xsf')
struc_40x40 = glob('/Users/nico/Desktop/simulation_outputs/percolation/40x40/structures/*xsf')

