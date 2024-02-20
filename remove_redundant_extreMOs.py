#!/usr/bin/env python

import sys
from os import path
import numpy as np
from utils_arpack import find_redundant_MOs


e_eps = 1e-4 #tolerance threshold for distinct energies

nn = sys.argv[1]
mo_types = ['lo', 'hi']

for mo_type in mo_types:

    e1 = np.load(path.expanduser(f'~/scratch/ArpackMAC/40x40/energies/{mo_type}/eARPACK_{mo_type}_bigMAC-{nn}.npy'))
    e2 = np.load(path.expanduser(f'~/scratch/ArpackMAC/40x40/energies/{mo_type}2/eARPACK_{mo_type}2_bigMAC-{nn}.npy'))
    M1 = np.load(path.expanduser(f'~/scratch/ArpackMAC/40x40/MOs/{mo_type}/MOs_ARPACK_{mo_type}_bigMAC-{nn}.npy'))
    M2 = np.load(path.expanduser(f'~/scratch/ArpackMAC/40x40/MOs/{mo_type}2/MOs_ARPACK_{mo_type}2_bigMAC-{nn}.npy'))

    rMOs = find_redundant_MOs(e1,e2,M1,M2)
    if rMOs is not None:
        print(rMOs)
        np.save(f'redundant_MO_inds_{mo_type}-{nn}.npy', np.vstack(rMOs))
    

