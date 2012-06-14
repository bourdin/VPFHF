#!/usr/bin/env python
import os

nx = 101
ny = 101
nz = 101
lx = 1.0
ly = 1.0
lz = 1.0
h = lx/nx
eta = 1e-8
#for k in (.25,.5,1,1.5,2,2.5,3):
#for k in (.5,1,1.5):
#for k in (2,2.5,3):
#for k in (2.5,3.0):
#for k in (.25,.25):

k = 1.5
eps = k * h
gc = 1./(1.+h/2./eps)
runcmd = 'echo qsub -l h_rt=01:00:00 -q development -pe 12way 240 -v EPS=%f,GC=%f,NX=%i,NY=%i,NZ=%iETA=%f ${VFDIR}/scripts/test16.sge'%(eps,gc,nx,ny,nz,eta)
os.system(runcmd)
