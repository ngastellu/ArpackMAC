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
        yield np.abs(np.diff(energies))

def plot_energy_spacing_stats(edir,lbls,filename_template='eARPACK_bigMAC',hist_kwargs=None):
    energies_ensemble = gen_energies(edir,lbls,filename_template)
    diffs = np.hstack(gen_energy_spacings(energies_ensemble))
    plt_utils.histogram(diffs,**hist_kwargs)


def find_redundant_MOs(e1,e2,M1,M2,e_eps=1e-4,M_eps=1e-5):
    """Given two sets of eigenpairs (e1, M1) and (e2, M2), return the eigenpairs that are common to both.
    This is done in two steps:

        1. Find the pairs of eigenvalues that within `e_eps` of (i.e. basically equal to) each other.

        2. Of those close eigenvalues, see if their corresponding eigenvectors colinear, or equivalently 
           if the abs value of their inner product is with M_eps of 1.
    
    Returns the indices of equivalent eigenpairs in both sets."""
    
    dE = e1[:,None] - e2
    ii, jj = (np.abs(dE) < e_eps).nonzero()
    if len(ii) > 0:
        M1i = M1[:,ii].T
        M2j = M2[:,jj]
        inner_products = np.abs((M1i @ M2j).diagonal())
        
        isame = (inner_products > 1-M_eps).nonzero()
        if len(isame) > 0:
            return ii[isame], jj[isame]
    return None

   
