#!/usr/bin/env python

import sys
from tb import TBsystem
from qcnico.coords_io import get_coords

rCC = 1.8

positions = get_coords('bigMAC_40x40-100_relaxed.xsf', dump=False)
tbsys = TBsystem(positions,rCC)
#tbsys.diagonalise()