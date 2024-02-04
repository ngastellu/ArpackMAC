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

occ_starts = np.zeros(lbls_virt.shape[0], 'int')
virt_starts = np.zeros(lbls_virt.shape[0], 'int')


for k, n in enumerate(lbls_virt):
    found_iHOMO = False
    found_iLUMO = False

    print(f'****** {n} ******') 
    odN = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/gap_check/occupied/odN-{n}.npy')
    vdN = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/gap_check/virtual/vdN-{n}.npy')

    odN = odN[odN > 0]
    vdN = vdN[vdN > 0]

    matching_ind = (all_lbls_virt == n).nonzero()[0]
    print('n = ', n)
    print('n\' = ', all_lbls_virt[matching_ind])
    N = Natoms[matching_ind]

    Nhalf = N // 2

    if odN.shape[0] > 0:
        odN -= (Nhalf-1)
        zero_inds = (odN == 0).nonzero()[0]
        if zero_inds.shape[0] == 0:
            print(odN)
        else:
            iHOMO = zero_inds[0]
            print('Found HOMO! --> iHOMO = ', iHOMO)
            found_iHOMO = True
            if iHOMO > 0:
                if odN[iHOMO-1] == 1:
                    iLUMO = iHOMO - 1
                    found_iLUMO = True
                    print('Found LUMO (from odN)! --> iLUMO = ', iLUMO)
                else:
                    print('Skipped LUMO: odN[iHOMO-1] = ', odN[iHOMO-1])
    
    if vdN.shape[0] > 0:
        vdN -= Nhalf
        zero_inds = (vdN == 0).nonzero()[0]
        if zero_inds.shape[0] == 0:
            print(vdN)
        else:
            iLUMO = zero_inds[0]
            found_iLUMO = True
            print('Found LUMO! --> iLUMO = ', iLUMO)
            if iLUMO > 1:
                if vdN[iLUMO-1] == -1:
                    iHOMO = iLUMO - 1
                    found_iHOMO = True
                    print('Found HOMO (from vdN)! --> iHOMO =', iHOMO)
                else:
                    'Skipped HOMO: vdN[iLUMO -1] = ', vdN[iLUMO-1]

    if found_iHOMO:
        nHOMO_found += 1        
    if found_iLUMO:
        nLUMO_found += 1
    if found_iHOMO and found_iLUMO:
        nboth_found += 1

    print('\n')

print('*********************')
print('Nb. HOMO found: = ', nHOMO_found)
print('Nb. LUMO found: = ', nLUMO_found)
print('Nb. both found: = ', nboth_found)