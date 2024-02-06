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

occ_starts = np.ones(lbls_virt.shape[0], 'int') * -1000
virt_starts = np.ones(lbls_virt.shape[0], 'int') * -1000



for k, n in enumerate(lbls_virt):

    print(f'****** {n} ******') 
    odN = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/gap_check/occupied/odN-{n}.npy')
    vdN = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/gap_check/virtual/vdN-{n}.npy')

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

    # Re-save separated occupied/virtual eigenpairs
    obtained_eoccs = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/eARPACK/occupied/eARPACK_bigMAC-{n}.npy')
    obtained_Moccs = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/MOs_ARPACK/occupied/MOs_ARPACK_bigMAC-{n}.npy')

    obtained_evirs = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/eARPACK/virtual/eARPACK_bigMAC-{n}.npy')
    obtained_Mvirs = np.load(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/MOs_ARPACK/virtual/MOs_ARPACK_bigMAC-{n}.npy')

    true_eoccs = obtained_eoccs[occ_starts[k]:] 
    true_Moccs = obtained_Moccs[:,occ_starts[k]:] 

    true_evirs = obtained_evirs[virt_starts[k]:]
    true_Mvirs = obtained_Mvirs[:,virt_starts[k]:]

    # Check if there any virtual states in my 'occupied' states that are absent from my virtual states
    # if occ_starts[k] > 0 and vdN[0] > 0:
    #     fake_occ_inds = np.arange(occ_starts[k]) # all states up to occ_start[k] are virtual
    #     odn = odN[fake_occ_inds] + 1
    #     missed_virt_inds = (odn < vdN[0]).nonzero()[0] # vdN[0] is always positive; all obtained virt states are actually virt 
    #     print(vdN[0])
    #     print(missed_virt_inds)

    #     if len(missed_virt_inds) > 0:
    #         print(f'Found {len(missed_virt_inds)} missed virtual states! vdN[0] = {vdN[0]}; odn = {odn}')

    #         imissed = np.argsort(obtained_eoccs[missed_virt_inds]) # eoccs are sorted in descending order, whereas evirs are sorted in ascending order; need to keepthings consistent

    #         missed_evirs = obtained_eoccs[imissed]
    #         missed_Mvirs = obtained_Moccs[:,imissed] 

    #         true_evirs = np.hstack((missed_evirs, obtained_evirs[virt_starts[k]:]))
    #         true_Mvirs = np.hstack((missed_Mvirs, obtained_Mvirs[:,virt_starts[k]:]))
    # else:
    #     true_evirs = obtained_evirs[virt_starts[k]:]
    #     true_Mvirs = obtained_Mvirs[:,virt_starts[k]:]


    np.save(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/eARPACK/occupied_cleaned/eARPACK_bigMAC-{n}.npy', true_eoccs)
    np.save(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/MOs_ARPACK/occupied_cleaned/MOs_ARPACK_bigMAC-{n}.npy', true_Moccs)

    np.save(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/eARPACK/virtual_cleaned/eARPACK_bigMAC-{n}.npy', true_evirs)
    np.save(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/MOs_ARPACK/virtual_cleaned/MOs_ARPACK_bigMAC-{n}.npy', true_Mvirs)



print('*********************')
print('Nb of "good" runs (i.e. even number of C atoms) = ',lbls_virt.shape[0])
print('Nb of structures for whom all obtained virtual MOs are actually virtual = ',(virt_starts == 0).sum())
print('Nb of structures for whom all obtained occupied MOs are actually occupied = ',(occ_starts == 0).sum())