#!/usr/bin/env python

import sys
import numpy as np

disp_file = sys.argv[-1]

D = np.loadtxt(disp_file)
work = 0.
for i in range(D.shape[0]-1):
    work += (D[i+1,0] - D[i,0]) * (D[i,2]+D[i+1,2]) / 2.

print("Work of pressure forces: %g"%(work * 1e-3))