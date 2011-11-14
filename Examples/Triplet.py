#!/usr/bin/env python

#PBS -l select=2:ncpus=4
#PBS -V
#PBS -m a
#PBS -j oe
#PBS -o Doublet
#PBS -q gmrs

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
		f.write('%s \t\t %s\n'%(key,str(D[key])))
	f.close()
  
def DictwriteJSON(D,filename):
	try:
		import json
	
		jsonfile = open(filename,'w')
		jsonfile.write(json.encoder.JSONEncoder().encode(D))
		jsonfile.flush()
		jsonfile.close()
	except ImportError:
		print 'JSON module not available, skipping DictJSONwrite'

def main():
	import os
	import os.path
	import time
	
	
	# Timestamp:
	print '###\n### Script started at %s\n###\n'%time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())
	Param = {}
	mygetenv(Param,'NX',10)
	mygetenv(Param,'NY',15)
	mygetenv(Param,'NZ',5)
	mygetenv(Param,'DX',1000/Param['NX'])
	mygetenv(Param,'DY',1500/Param['NY'])
	mygetenv(Param,'DZ',100/Param['NZ'])
	mygetenv(Param,'TEMPR1',180)
	mygetenv(Param,'WTEMPR',80)
	mygetenv(Param,'PREFIX','Triplet')

	### PBS stuff
	mygetenv(Param,'PBS_JOBID','00000')
	mygetenv(Param,'PBS_JOBNAME','JOBNAME')
	mygetenv(Param,'PBS_O_WORKDIR',os.getenv('PWD'))

	mygetenv(Param,'VFDIR')
	mygetenv(Param,'GMRSDIR')
	mygetenv(Param,'GMRSARCH')
	mygetenv(Param,'PETSC_DIR')
	mygetenv(Param,'PETSC_ARCH')
	mygetenv(Param,'GMRSBIN','GMRS_VF')
	
	mygetenv(Param,'MODE','ELASTICITY')
	mygetenv(Param,'PRESET','SYMXY')
	mygetenv(Param,'COUPLING','FULL')
	mygetenv(Param,'VERBOSE',0)

	mygetenv(Param,'E',1e+6)
 	mygetenv(Param,'EPSILON',20)
	mygetenv(Param,'NU',.25)
	mygetenv(Param,'GC',1)
	mygetenv(Param,'ALPHA',1e-5)
	mygetenv(Param,'ETA',1e-5)
	mygetenv(Param,'INSITUMIN','0,0,0')
	mygetenv(Param,'INSITUMAX','0,0,0')
	mygetenv(Param,'VFOPTS','-U_ksp_max_it 20000 -V_ksp_max_it 2000')

	print 'Param:    \n',Param

	workdir = os.path.join(Param['PBS_O_WORKDIR'],'Doublet',Param['PBS_JOBID'])
	print 'Work dir is %s'%workdir
	if not os.path.exists(workdir):
		os.makedirs(workdir)
	os.chdir(workdir)

	### Save computation parameters to a txt and a json file
	Dictwritetxt(Param,'00_INFO.txt')
	DictwriteJSON(Param,'00_INFO.json')
	
	### Generate GMRS data file
	filenamein  = os.path.join(Param['PBS_O_WORKDIR'],Param['PREFIX']+'.dat')
	filenameout = '%s.dat'%Param['PBS_JOBID']
	print filenamein,filenameout
	open(filenameout,'w').write( open(filenamein,'r').read()%Param ) 
	
	open('temp.txt','w').write('%s\n%s.out\n'%(filenameout,Param['PBS_JOBID']))

	###
	### Run the computation
	###
	
	t1 = time.time()
	cmd = '''mpirun %(GMRSDIR)s/%(GMRSARCH)s/%(GMRSBIN)s -p %(PREFIX)s \
           -E %(E)f -nu %(NU)f -alpha %(ALPHA)f -gc %(GC)f -mode %(MODE)s \
           -preset %(PRESET)s -coupling %(COUPLING)s -epsilon %(EPSILON)f -eta %(ETA)f \
           -verbose %(VERBOSE)i -insitumin %(INSITUMIN)s -insitumax %(INSITUMAX)s %(VFOPTS)s < temp.txt'''%Param
	print cmd
	os.system(cmd)
	t2 = time.time()
	print '###\n### Computation took %fs\n###\n'%(t2-t1)
	###
	### Convert petsc binary files to hdf5 format
	### This needs to be done separately as parallel hdf5 
	### does not work with scali mpi and pgi compilers
	###
	cmd = 'mpirun -np 1 %(VFDIR)s/bin/%(PETSC_ARCH)s/h5export -p %(PREFIX)s'%Param
	print cmd
	#os.system(cmd)
	print '###\n### Script finished at %s\n###\n'%time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

import sys  
if __name__ == "__main__":
	sys.exit(main())

