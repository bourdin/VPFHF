#!/usr/bin/env python
import os

nx = 355
ny = 355
lx = 8.0
h = lx/nx
eta = 1e-8
theta=20
#for k in (.25,.5,1,1.5,2,2.5,3):
#for k in (.5,1,1.5):
#for k in (2,2.5,3):
for k in (1.,1.5,2.):
#for k in (.25,.25):
    eps = k * h
    gc = 1./(1.+h/2./eps)
    runcmd = 'qsub -l h_rt=01:00:00 -q development -pe 12way 192 -v EPS=%f,GC=%f,NX=%i,NY=%i,ETA=%f,THETA=%f ${VFDIR}/scripts/test28.sge'%(eps,gc,nx,ny,eta,theta)
    os.system(runcmd)
