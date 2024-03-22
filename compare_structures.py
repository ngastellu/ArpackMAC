#!/usr/bin/env python

from collections import deque
from glob import glob
import numpy as np
from scipy.spatial import KDTree
from qcnico.coords_io import read_xyz, read_xsf
from qcnico.remove_dangling_carbons import remove_dangling_carbons
from qcnico.plt_utils import multiple_histograms



def get_nn_dists(posfile):
    rCC = 1.8
    if posfile.split('.')[-1] == 'xyz':
        pos = remove_dangling_carbons(read_xyz(posfile),rCC)
    elif posfile.split('.')[-1] == 'xsf':
        pos = remove_dangling_carbons(read_xsf(posfile)[0],rCC)
    tree = KDTree(pos)
    M = tree.sparse_distance_matrix(tree,rCC)
    dists = np.array(list(M.values()))
    dists = np.unique(dists[dists>0])
    return dists



strucdir = '/Users/nico/Desktop/simulation_outputs/MAC_structures/'

t1_unrel_strucs = glob(strucdir + 'Ata_structures/t1/*xyz')
t1_rel_strucs = glob(strucdir + 'Ata_structures/t1/relaxed/*xsf')
struc_40x40 = glob('/Users/nico/Desktop/simulation_outputs/percolation/40x40/structures/*xsf')


d40x40 = deque()
drel = deque()

print('Building 40x40 deque...')
for f in struc_40x40:
    d40x40.append(get_nn_dists(f))
print('Done!')

print('Building relaxed deque...')
for f in t1_rel_strucs:
    drel.append(get_nn_dists(f))
print('Done!')


dists_40x40 = np.hstack(d40x40)
dists_rel = np.hstack(drel)

multiple_histograms((dists_40x40,dists_rel), ('Michael-MAC', 'AtaT1-MAC'), nbins=1000, xlabel='Bond length [$\\text{\AA}$]', normalised=True)