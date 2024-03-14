#!/usr/bin/env python

import numpy as np


hdir_relaxed = '/Users/nico/Desktop/simulation_outputs/percolation/Ata_structures/t1/Hao_ARPACK/'
hdir_unrelaxed = '/Users/nico/Desktop/simulation_outputs/percolation/Ata_structures/t1/Hao_ARPACK_unrelaxed/'

labels = set(range(38)) - {4,16,18,21,27,32} #remove structures for which the diagonlisation was successful when unrelaxed

for n in labels:
    ii_relaxed = np.load(hdir_relaxed + f'inds/ii-{n}.npy')
    jj_relaxed = np.load(hdir_relaxed + f'inds/jj-{n}.npy')

    ii_unrelaxed = np.load(hdir_unrelaxed + f'inds/ii-{n}.npy')
    jj_unrelaxed = np.load(hdir_unrelaxed + f'inds/jj-{n}.npy')

    hvals_relaxed = np.load(hdir_relaxed + f'hvals/hvals-{n}.npy')
    hvals_unrelaxed = np.load(hdir_unrelaxed + f'hvals/hvals-{n}.npy')

    print(f'\n***** {n} *****')
    print('ii match = ', np.all(ii_relaxed == ii_unrelaxed))
    print('jj match = ', np.all(jj_relaxed == jj_unrelaxed))
    print('hvals match = ', np.all(hvals_relaxed == hvals_unrelaxed))
    
    print("sum{ii_relaxed - ii_unrelaxed} = ", np.sum(ii_relaxed - ii_unrelaxed))
    print("sum{jj_relaxed - jj_unrelaxed} = ", np.sum(jj_relaxed - jj_unrelaxed))
    print("sum{hvals_relaxed - hvals_unrelaxed} = ", np.sum(hvals_relaxed - hvals_unrelaxed))

