#!/usr/bin/env python

import numpy as np
import os


def check_dE(energies, dE_max):
    energies = np.sort(energies)
    eHOMO = energies[0]
    eLUMO = energies[1]
    eF = 0.5 * (eHOMO + eLUMO)
    return np.abs(energies[1:] - eF) < dE_max

kB = 8.617333e-5
T = 300
prefactor = 4.0
dE_max = prefactor * kB * T


edir = 'energies/virtual_w_HOMO'

for efile in os.listdir(edir):
    n = int(efile.split('-')[-1].split('.')[0])
    print(f'{n} --> ', end='')
    energies = np.load(os.path.join(edir, efile))
    check = check_dE(energies, dE_max)
    if not np.all(check):
        print(f'{(~check).sum()} / {check.shape[0]} problematic energies')
    else:
        print('all g!')
        
