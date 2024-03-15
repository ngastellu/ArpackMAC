#!/usr/bin/env python

import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from glob import glob
from qcnico.plt_utils import setup_tex, multiple_histograms, histogram

nstrucs = 38

Hdir_relaxed = '/Users/nico/Desktop/simulation_outputs/percolation/Ata_structures/t1/Hao_ARPACK/hvals/'
Hdir_unrelaxed = '/Users/nico/Desktop/simulation_outputs/percolation/Ata_structures/t1/Hao_ARPACK_unrelaxed/hvals/'

hvals_relaxed = np.hstack([np.load(Hdir_relaxed + f'hvals-{n}.npy') for n in range(nstrucs)])
hvals_unrelaxed = np.hstack([np.load(Hdir_unrelaxed + f'hvals-{n}.npy') for n in range(nstrucs)])

    
multiple_histograms((hvals_relaxed,hvals_unrelaxed),('relaxed','unrelaxed'),nbins=500, xlabel='$H_{ij}$ [eV]', log_counts=True)

histogram(hvals_relaxed,nbins=1000, xlabel='$H_{ij}$ [eV]')

del hvals_unrelaxed

Hdir_40x40 = '/Users/nico/Desktop/simulation_outputs/percolation/40x40/Hao_ARPACK/hvals/'
hvals_40x40 = np.hstack([np.load(f)] for f in glob(Hdir_40x40 + '*.npy')) 
multiple_histograms((hvals_40x40, hvals_relaxed),('Micheal MAC', 'AtaT1 MAC'),nbins=1000,xlabel='$H_{ij}$ [eV]')

