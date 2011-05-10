#!/usr/bin/env python

#PBS -l nodes=2:ppn=2
#PBS -V
#PBS -m ae
#PBS -M bnvk@chevron.com
#PBS -j oe

def mygetenv(Dict,key,defaultvalue=None):
  ### I could do better 
  import os
  tmp = os.getenv(key)
  if tmp == None:
    Dict[key] = defaultvalue
  else:
    try:
      Dict[key] = float(tmp)
    except ValueError:
      Dict[key] = tmp

def Dictwritetxt(D,filename):
  f = open(filename,'a')
  K = D.keys()
  K.sort()
  for key in K:
    f.write('%s \t\t %s\n'%(key, str(D[key])))
  f.close()
  


def main():
	import os
	import os.path
	import numpy as np
	import hashlib
	import shutil
	import json
	
	Param = {}
	mygetenv(Param,'NX',10)
	mygetenv(Param,'NY',10)
	mygetenv(Param,'NZ',5)
	mygetenv(Param,'DX',100)
	mygetenv(Param,'DY',100)
	mygetenv(Param,'DZ',20)
	mygetenv(Param,'TEMPR1',180)
	mygetenv(Param,'WTEMPR',80)
	mygetenv(Param,'PREFIX','Doublet')

	### PBS stuff
	mygetenv(Param,'PBS_JOBID','00000')
	mygetenv(Param,'PBS_JOBNAME','JOBNAME')
	mygetenv(Param,'PBS_O_WORKDIR',os.getenv('PWD'))
	mygetenv(Param,'MODE','GMRS_FakeVF')

	mygetenv(Param,'VFDIR')
	mygetenv(Param,'GMRSDIR')
	mygetenv(Param,'GMRSARCH')
	mygetenv(Param,'GMRSBIN')
	 
	
	print 'Param:    \n', Param

	workdir = os.path.join(Param['PBS_O_WORKDIR'],Param['PREFIX']+'-'+Param['PBS_JOBID'])
	print 'Work dir is %s'%workdir
	if not os.path.exists(workdir):
		os.makedirs(workdir)
	os.chdir(workdir)

	### Save computation parameters to a txt and a json fil
	Dictwritetxt(Param,'00_INFO.txt')
	
	jsonfile = open('00_INFO.json','w')
	jsonfile.write(json.encoder.JSONEncoder().encode(Param))
	jsonfile.flush()
	jsonfile.close()
	
	### Generate GMRS data file
	filenamein  = os.path.join(Param['PBS_O_WORKDIR'],Param['PREFIX']+'.dat')
	filenameout = '%s.dat'%Param['PBS_JOBID']
	print filenamein, filenameout
	open(filenameout,'w').write( open(filenamein,'r').read()%Param ) 
	
	open('temp.txt','w').write('%s\n%s.out\n'%(filenameout,Param['PBS_JOBID']))
	
	cmd = 'mpirun %s < temp.txt'%os.path.join(Param['GMRSDIR'],Param['GMRSARCH'],Param['GMRSBIN'])
	print cmd
	#os.system(cmd)
	
import sys  
if __name__ == "__main__":
	sys.exit(main())

