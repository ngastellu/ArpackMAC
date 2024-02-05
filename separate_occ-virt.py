#!/usr/bin/env python 

import numpy as np
from glob import glob

occfiles = glob('/Users/nico/Desktop/simulation_outputs/percolation/40x40/gap_check/occupied/*npy')
lbls_occ = np.sort([int(f.split('-')[1].split('.')[0]) for f in occfiles])

virtfiles = glob('/Users/nico/Desktop/simulation_outputs/percolation/40x40/gap_check/occupied/*npy')
lbls_virt = np.sort([int(f.split('-')[1].split('.')[0]) for f in virtfiles])

all_lbls_virt = np.load('vir_lbls.npy')

print(np.all(lbls_occ == lbls_virt)) #True

Natoms = np.load('natoms.npy')

nHOMO_found = 0
nLUMO_found = 0
nboth_found = 0

occ_starts = np.ones(lbls_virt.shape[0], 'int') * -1000
virt_starts = np.ones(lbls_virt.shape[0], 'int') * -1000



for k, n in enumerate(lbls_virt):
    found_iHOMO = False
    found_iLUMO = False

    print(f'****** {n} ******') 
    odN = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/gap_check/occupied/odN-{n}.npy')
    vdN = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/gap_check/virtual/vdN-{n}.npy')

    odN = odN
    vdN = vdN

    matching_ind = (all_lbls_virt == n).nonzero()[0]
    N = Natoms[matching_ind]

    Nhalf = N // 2

    odN -= (Nhalf-1)
    occ_inds = (odN <= 0).nonzero()[0]
    occ_starts[k] = np.min(occ_inds)
    print(occ_starts[k])

    vdN -= Nhalf
    virt_inds = (odN >= 0).nonzero()[0]
    virt_starts[k] = np.min(virt_inds)
    print(virt_starts[k])

    print('\n')

print('*********************')
print('Nb of "good" runs (i.e. even number of C atoms) = ',lbls_virt.shape[0])
print('Nb of structures for whom all obtained virtual MOs are actually virtual = ',(virt_starts == 0).sum())
print('Nb of structures for whom all obtained occupied MOs are actually occupied = ',(occ_starts == 0).sum())