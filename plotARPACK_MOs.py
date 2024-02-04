from os import path
import numpy as np
import matplotlib.pyplot as plt
from qcnico.coords_io import read_xsf
from qcnico.qcplots import plot_MO


datadir = path.expanduser("~/Desktop/simulation_outputs/percolation/40x40")
posdir = path.join(datadir, "structures")

lbls = [f'bigMAC-{n}' for n in [2,3,5,6,7,9,10]]

def getARPACKdata(datadir,lbls):
    edir = path.join(datadir,'eARPACK')
    Mdir = path.join(datadir,'MOs_ARPACK')
    ee = [np.load(path.join(edir,f'eARPACK_{lbl}.npy')) for lbl in lbls]
    MM = [np.load(path.join(Mdir,f'MOs_ARPACK_{lbl}.npy')) for lbl in lbls]
    return ee, MM

def get_pos(posdir,lbls):
    return [read_xsf(path.join(posdir,f"{lbl}_relaxed.xsf"))[0] for lbl in lbls]
    


ee, MM = getARPACKdata(datadir,lbls)
posarrs = get_pos(posdir,lbls)

for pos, M in zip(posarrs,MM):
    plot_MO(pos,M,0,dotsize=1.0,show_COM=True,show_rgyr=True)