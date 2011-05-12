#!/usr/bin/env python

#PBS -l nodes=2:ppn=2
#PBS -V
#PBS -m ae
#PBS -M bnvk@chevron.com
#PBS -j oe
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

	mygetenv(Param,'VFDIR')
	mygetenv(Param,'GMRSDIR')
	mygetenv(Param,'GMRSARCH')
	mygetenv(Param,'PETSC_DIR')
	mygetenv(Param,'PETSC_ARCH')
	mygetenv(Param,'GMRSBIN','GMRS_VF')
	
	mygetenv(Param,'MODE','ELASTICITY')
	mygetenv(Param,'PRESET','SIMXY')

	mygetenv(Param,'E',5e+3)
	mygetenv(Param,'NU',.25)
	mygetenv(Param,'GC',5e-2)
	mygetenv(Param,'ALPHA',1e-5)
	mygetenv(Param,'VFOPTS','-U_ksp_monitor -V_ksp_monitor -verbose 1 -U_ksp_view -V_ksp_view')

	print 'Param:    \n',Param

	workdir = os.path.join(Param['PBS_O_WORKDIR'],Param['PREFIX']+'-'+Param['PBS_JOBID'])
	print 'Work dir is %s'%workdir
	if not os.path.exists(workdir):
		os.makedirs(workdir)
	os.chdir(workdir)

	### Save computation parameters to a txt and a json file
	Dictwritetxt(Param,'00_INFO.txt')
	DictwriteJSON(Param,'0_INFO.json')
	
	### Generate GMRS data file
	filenamein  = os.path.join(Param['PBS_O_WORKDIR'],Param['PREFIX']+'.dat')
	filenameout = '%s.dat'%Param['PBS_JOBID']
	print filenamein,filenameout
	open(filenameout,'w').write( open(filenamein,'r').read()%Param ) 
	
	open('temp.txt','w').write('%s\n%s.out\n'%(filenameout,Param['PBS_JOBID']))

	###
	### Run the computation
	###
	
	cmd = 'mpirun %(GMRSDIR)s/%(GMRSARCH)s/%(GMRSBIN)s -p %(PREFIX)s -E %(E)f -nu %(NU)f -alpha %(ALPHA)f -gc %(GC)f -mode %(MODE)s -preset %(PRESET)s %(VFOPTS)s < temp.txt'%Param
	print cmd
	os.system(cmd)
	
	###
	### Convert petsc binary files to hdf5 format
	### This needs to be done separately as parallel hdf5 
	### does not work with scali mpi and pgi compilers
	###
	cmd = 'mpirun -np 1 %s -p %s'%(os.path.join(Param['VFDIR'],'bin','h5export'),Param['PREFIX'])
	print cmd
	os.system(cmd)
	
import sys  
if __name__ == "__main__":
	sys.exit(main())

