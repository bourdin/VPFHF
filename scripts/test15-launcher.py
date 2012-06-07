#!/usr/bin/env python
import os

nx = 1000.0
lx = 8.0
h = lx/nx
for k in (.5,1,1.5,2,2.5,3):
    eps = k * h
    gc = 1./(1.+h/4./eps)
    runcmd = 'echo qsub -v EPS=%f,GC=%f ./HS2D.sge'%(eps,gc)
    print runcmd
    os.system(runcmd)
