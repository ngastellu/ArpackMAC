#!/usr/bin/env python

from os import path
import numpy as np
import matplotlib.pyplot as plt
from qcnico.coords_io import read_xyz
from qcnico.qcffpi_io import load_qcffpi_data
from qcnico.plt_utils import setup_tex, histogram
from utils_arpackMAC import find_redundant_MOs

"""This script compares the results of ArpackMAC with those of QCFFPI-TB for a single MAC structure containing ~6000 atoms.
QCFFPI-TB obtains the full eigenspectrum of a given Hamiltonian while ArpackMAC, only obtains the lowest- and highest-energy MOs,
along with the low-lying virtual MOs."""



posfile = '/Users/nico/Desktop/simulation_outputs/MAC_structures/truncated_40x40-101.xyz'
pos = read_xyz(posfile)


datadir = '/Users/nico/Desktop/simulation_outputs/ArpackMAC_vs_qcffpiTB'

# Load QCFFPI-TB data
qcffpi_datadir = path.join(datadir,"qcffpiTB")
orbfile = path.join(qcffpi_datadir,"orb_energy.dat")
Mfile_q = path.join(qcffpi_datadir,"MO_coefs.dat") 
pos_q, energies_q, M_q = load_qcffpi_data(orbfile,Mfile_q)

print('Position arrays match: ', np.all(pos == pos_q))

# Load ArpackMAC data
motypes = ['lo', 'virtual_w_HOMO', 'hi']
arpackMAC_datadir = path.join(datadir,"ArpackMAC")
energies_a = np.hstack([np.load(path.join(arpackMAC_datadir,f"eARPACK_{mt}_40x40-101.npy")) for mt in motypes])
M_a = np.hstack([np.load(path.join(arpackMAC_datadir,f"MOs_ARPACK_{mt}_40x40-101.npy")) for mt in motypes])


imatch = find_redundant_MOs(energies_q, energies_a, M_q, M_a)

print(imatch)

ematch_q = energies_q[imatch[0]]
Mmatch_q = M_q[imatch[0]]

ematch_a = energies_a[imatch[1]]
Mmatch_a = M_a[imatch[1]]

# for y = x line
glob_min = np.min(np.hstack(ematch_q,ematch_a)) - 0.01
glob_max = np.max(np.hstack(ematch_q,ematch_a)) + 0.01
x = np.linspace(glob_min, glob_max, 100)

plt.plot(ematch_q, ematch_q, 'ro')
plt.plot(x,x,'k--') # y = x line to guide the eye
plt.xlabel('Exact energies [eV]')
plt.ylabel('Approximate energies [eV]')
plt.show()


ediffs = ematch_q - ematch_a
histogram(ediffs,bins=50,xlabel='$E_{\\text{exact}} - E_{\\text{approx}}$')

dotprods = np.diagonal((Mmatch_q.T) @ Mmatch_a.T)
histogram(ediffs,bins=50,xlabel='$\langle\psi_{\\text{exact}}|\psi_{\\text{approx}}\\rangle$')