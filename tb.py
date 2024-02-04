#!/usr/bin/env python

import numpy as np
from qcnico.graph_tools import adjacency_matrix_sparse
from remove_dangling_carbons import remove_dangling_carbons



def Htb(pos, rNN, t=1):
    _, neighbour_pairs = adjacency_matrix_sparse(pos,rNN,return_pairs=True)
    N = pos.shape[0]
    H = np.zeros((N,N),dtype=type(t))
    rows, cols = neighbour_pairs.T
    H[rows,cols] = t
    H += H.T #symmetrise Hamiltonian

    return H

def Htb_lindberg(pos, rNN):
    '''Constructs tight-binding Hamiltonian in which the hopping parameters are distance-dependent, and computed using the Lindberg approximation.
    The parameter values are obtained from table III in Warshel and Karplus, JACS, 94, 5612 (1972).'''

    #Parameters for resonance integral between two conjugated C atoms
    b0 = -2.438 #eV
    kb = 0.405 #angstrom^-1
    R0 = 1.397 #angstrom
    mub = 2.035 #angstrom^-1

    N = pos.shape[0]
    dist_mat = np.linalg.norm(pos[None,:] - pos[:,None], axis=2)
    print('dist_mat.shape = ', dist_mat.shape)
    np.fill_diagonal(dist_mat,100) # do this to ignore diagonal entries of distance matrix (trivially zero)
    inds = (dist_mat < rNN).nonzero()
    dnns = dist_mat[dist_mat < rNN]
    print('dnns.shape = ', dnns.shape)
    dnns -= R0
    dnns = b0*np.exp(-mub*dnns)*(1 + kb*dnns)
    H = np.zeros((N,N),dtype=float)
    H[inds] = dnns
    return H


class TBsystem:
    def __init__(self, pos, rNN, t=1, remove_dangling=True, lindberg=False):
        self.rNN = rNN
        self.t = t
        if remove_dangling:
            self.pos = remove_dangling_carbons(pos, rNN)
        else:
            self.pos = pos
        if lindberg:
            self.H = Htb_lindberg(self.pos,rNN)
        else:
            self.H = Htb(self.pos, rNN, t)

    def diagonalise(self):
        ee, M = np.linalg.eig(self.H)
        sorted_inds = np.argsort(ee)

        # sort eigenvalues and eigenvectors in order of increasing eigenvalues
        self.energies = ee[sorted_inds]
        self.M = M[:,sorted_inds]
