#!/usr/bin/env python

#PBS -l nodes=2:ppn=2
#PBS -V
#PBS -m a
#PBS -M bnvk@chevron.com
#PBS -j oe
#PBS -o FakeGMRS_VF

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
    #mygetenv(Param,'NX',10)
    #mygetenv(Param,'NY',10)
    #mygetenv(Param,'NZ',10)
    #mygetenv(Param,'DX',1000/Param['NX'])
    #mygetenv(Param,'DY',1000/Param['NX'])
    #mygetenv(Param,'DZ',200/Param['NX'])
    #mygetenv(Param,'TEMPR1',180)
    #mygetenv(Param,'WTEMPR',80)
    #mygetenv(Param,'PREFIX','Doublet')

    ### PBS stuff
    mygetenv(Param,'PBS_JOBID','00000')
    mygetenv(Param,'PBS_JOBNAME','JOBNAME')
    mygetenv(Param,'PBS_O_WORKDIR',os.getenv('PWD'))
    mygetenv(Param,'PETSC_ARCH',os.getenv('PETSC_ARCH'))
    mygetenv(Param,'VFDIR',os.getenv('VFDIR'))
    mygetenv(Param,'MODE','FakeGMRS_VF')
<<<<<<< /Volumes/VTrak1/Users/bourdin/Development/VF-Chevron/Examples/FakeGMRS_VF.py
    mygetenv(Param,'U','-U_ksp_max_it 20000')
    mygetenv(Param,'V','-V_ksp_max_it 20000')
=======
    mygetenv(Param,'U',' ')
    mygetenv(Param,'V','')
>>>>>>> /tmp/FakeGMRS_VF.py~other.AL6xwP
    mygetenv(Param,'FMT','hdf5')
    mygetenv(Param,'OPTS','-alpha 1e-2 -gc .1 -epsilon 20 -eta 1e-5 -mode FRACTURE -preset SYMXY')

    #mygetenv(Param,'FAKEGMRSBIN',os.path.join(os.getenv('VFDIR'),os.getenv('PETSC_ARCH'),os.getenv('MODE')))
    
    print 'Param:    \n',Param

    workdir = os.path.join(Param['PBS_O_WORKDIR'],'FakeGMRS_VF',Param['PBS_JOBID'])
    print 'Work dir is %s'%workdir
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    os.chdir(workdir)

    ### Save computation parameters to a txt and a json file
    Dictwritetxt(Param,'00_INFO.txt')
    DictwriteJSON(Param,'00_INFO.json')
    
    ###
    ### Run the computation
    ###
    
    t1 = time.time()
    cmd = 'mpirun %(VFDIR)s/bin/%(PETSC_ARCH)s/%(MODE)s %(U)s %(V)s -format %(FMT)s %(OPTS)s -p %(PBS_JOBID)s'%Param
    print cmd
    os.system(cmd)
    t2 = time.time()
    print '###\n### Computation took %fs\n###\n'%(t2-t1)
    print '###\n### Script finished at %s\n###\n'%time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())

import sys  
if __name__ == "__main__":
    sys.exit(main())

