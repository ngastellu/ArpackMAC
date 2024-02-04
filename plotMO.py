#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from qcnico.qcplots import plot_MO
from qcnico import plt_utils


M = np.load("MOs_ARPACK_bigMAC_iLUMO=1.npy")
pos  = np.load('pos.npy')

pos[:,0] -= np.min(pos[:,0])
pos[:,1] -= np.min(pos[:,1])

plt_utils.setup_tex()

for i in range(4):
    plot_MO(pos, M, i)
    