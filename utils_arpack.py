#!/usr/bin/env python

import numpy as np
from qcnico.coords_io import read_xsf
from qcnico.remove_dangling_carbons import remove_dangling_carbons
from qcnico import plt_utils


def gen_mos(Mdir,lbls,filename_template='MOs_ARPACK_bigMAC'):
    for nn in lbls:
        yield np.load(Mdir+f'/{filename_template}-{nn}.npy')

def gen_energies(edir,lbls,filename_template='eARPACK_bigMAC'):
    for nn in lbls:
        yield np.load(edir+f'/{filename_template}-{nn}.npy')


def gen_pos(posdir, lbls, rCC):
    for nn in lbls:
        yield remove_dangling_carbons(read_xsf(path.join(posdir,f"bigMAC-{nn}_relaxed.xsf"))[0], rCC)

def gen_energy_spacings(energies_ensemble):
    for energies in energies_ensemble:
        yield np.diff(energies_ensemble)

def plot_energy_spacing_stats(edir,lbls,filename_template='eARPACK_ARPACK_bigMAC',hist_kwargs=None):
    energies_ensemble = gen_energies(edir,lbls,filename_template)
    diffs = np.hstack(gen_energy_spacings(energies_ensemble))
    plt_utils.histogram(diffs,**hist_kwargs)
   
