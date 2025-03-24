#!/usr/bin/env python

import sys
from os import path
import numpy as np
from utils_arpack import find_redundant_MOs
from qcnico.qcplots import plot_MO
from qcnico.coords_io import read_xsf


e_eps = 1e-4 #tolerance threshold for distinct energies

# nn = sys.argv[1]
nn = 10
mo_type = 'lo'

posdir = '/Users/nico/Desktop/simulation_outputs/percolation/40x40/structures/'

# for mo_type in mo_types:

e1 = np.load(path.expanduser(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/eARPACK/{mo_type}/eARPACK_{mo_type}_bigMAC-{nn}.npy'))
e2 = np.load(path.expanduser(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/eARPACK/{mo_type}2/eARPACK_{mo_type}2_bigMAC-{nn}.npy'))
M1 = np.load(path.expanduser(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/MOs_ARPACK/{mo_type}/MOs_ARPACK_{mo_type}_bigMAC-{nn}.npy'))
M2 = np.load(path.expanduser(f'/Users/nico/Desktop/simulation_outputs/percolation/40x40/MOs_ARPACK/{mo_type}2/MOs_ARPACK_{mo_type}2_bigMAC-{nn}.npy'))

rMOs = find_redundant_MOs(e1,e2,M1,M2)
if rMOs is not None:
    print(rMOs)
    rMOs = np.vstack(rMOs).T
    print(rMOs.shape)
    # pos, _ = read_xsf(posdir+f'bigMAC-{nn}_relaxed.xsf')
    # for ij in rMOs:
    #     i, j = ij
    #     print(f'**** ({i},{j}) ****')
    #     print(f'Energy difference = {np.abs(e1[i] - e2[j])} eV')
    #     print(f'|1 - <i|j> | = {np.abs(1 - np.dot(M1[:,i],M2[:,j]))}')
    #     # plot_MO(pos,M1,i,dotsize=0.1,cmap='coolwarm', plot_amplitude=True, show_rgyr=True)
    #     # plot_MO(pos,M2,j,dotsize=0.1,cmap='coolwarm', plot_amplitude=True, show_rgyr=True)



