#!/usr/bin/env python

def parseCommandLine(debug=False):
    import argparse 
    import random
    parser = argparse.ArgumentParser()
    ### Mesh stuff
    vfArgs = parser.add_argument_group('VF_Standalone')
    vfArgs.add_argument('--n',type=int,nargs=3,default=[100,2,100],help='number of grid points')
    vfArgs.add_argument('--l',type=float,nargs=3,default=[8,.1,8],help='geometry')
    vfArgs.add_argument('--altmintol',type=float,default=1e-4,help='Tolerance on alternate minimization')
    vfArgs.add_argument('--altminmaxit',type=int,default=100,help='Maximum number of alternate minimizations')
    vfArgs.add_argument('--atnum',type=int,default=1,choices=[1,2],help='AT variant')
    vfArgs.add_argument('--epsilon',type=float,default=None,help='AT regularization length')
    vfArgs.add_argument('--eta',type=float,default=None,help='AT residual stiffness')
    vfArgs.add_argument('--format',default='vtk',choices=['bin','hdf5','vtk'],help='file format')

    vfArgs.add_argument('--E',type=float,default=1.,help='Young modulus')
    vfArgs.add_argument('--nu',type=float,default=0.,help='Poisson ratio')
    vfArgs.add_argument('--gc',type=float,default=1.,help='Fracture toughness')
    vfArgs.add_argument('--gceff',default=False,action='store_true',help='Modify Gc to account for effective toughness')
    vfArgs.add_argument('--prefix',help='Name files after prefix + geometry instead of job ID',default=None)

    vfArgs.add_argument('--insitumin',type=float,nargs=6,default=[0,0,0,0,0,0],help='Insitu stresses in the plane Z0')
    vfArgs.add_argument('--insitumax',type=float,nargs=6,default=None,help='Insitu stresses in the plane Z1')
    vfArgs.add_argument('--gamg',default=False,action='store_true',help='Use PETSc GAMG solvers')
    
    misc = parser.add_argument_group('misc')
    misc.add_argument('--mpiexec',help='mpi exec command',default='srun')
    misc.add_argument('--extraopts',help='additional options passed to test16',default='')
    test16 = parser.add_argument_group('test16')
    test16.add_argument('--minvol',type=float,default=0.,help='Initial injected volume')
    test16.add_argument('--maxvol',type=float,default=.1,help='Final injected volume')
    test16.add_argument('--maxtimestep',type=int,default=25,help='Number of time steps')    
    test16.add_argument('--mode',type=int,default=3,help='BC type')
    
    
    pc = parser.add_argument_group('penny-shaped cracks')
    pc.add_argument('--npc',type=int,default=1,help='Number of penny shaped cracks')

    args,CL = parser.parse_known_args()
    
    h = min( [l/float(n) for l,n in zip(args.l,args.n)])
    for c in range(args.npc):
      pc.add_argument('--pc%i_R'%c,type=float,default=.5,help='Radius of penny shaped cracks %i'%c)  
      pc.add_argument('--pc%i_center'%c,type=float,nargs=3,default=[x/2. for x in args.l],help='Center of penny shaped cracks %i'%c)  
      pc.add_argument('--pc%i_phi'%c,type=float,default=45,help='Colatitude of penny shaped cracks %i'%c)  
      pc.add_argument('--pc%i_theta'%c,type=float,default=0,help='Polar angle of penny shaped cracks %i'%c)  
      pc.add_argument('--pc%i_thickness'%c,type=float,default=1.5*h,help='Thickness of penny shaped cracks %i'%c)  
    

    args = parser.parse_args()
    if not args.epsilon:
      args.epsilon = 2*h
    if not args.eta:
      args.eta = 1.e-8 * args.epsilon
      
    if args.gceff:
      if args.atnum == 1:
        args.gc = args.gc * 1./(1.+3.*h/8./args.epsilon)
      elif args.atnum == 2:
        args.gc = args.gc * 1./(1.+h/2./args.epsilon)

    if not args.insitumax:
      args.insitumax = args.insitumin  

    args = vars(args)    
    return args
    
def main():
  import pprint
  import os
  import os.path
  import shutil
  import json
  
  
  args = parseCommandLine() 
  
  # Get prefix
  if os.getenv('PBS_JOBID'):
      args['prefix'] = os.getenv('PBS_JOBID')
  elif os.getenv('JOB_ID'):
      args['prefix'] = os.getenv('JOB_ID')
  elif os.getenv('SLURM_JOB_ID'):
      args['prefix'] = os.getenv('SLURM_JOB_ID')
  elif not args['prefix']:
      args['prefix']='test16'
  
  # Get workdir
  if os.getenv('PBS_O_WORKDIR'):
      # We are runnning inside a PBS job 
      args['workdir'] = os.path.join(os.getenv('PBS_O_WORKDIR'),args['prefix'])
  elif os.getenv('SGE_O_WORKDIR'):
      # We are running inside a SGE job
      args['workdir'] = os.path.join(os.getenv('SGE_O_WORKDIR'),args['prefix'])
  elif os.getenv('SLURM_SUBMIT_DIR'):
      # We are running inside a SBATCH / SRUN job
      args['workdir'] = os.path.join(os.getenv('SLURM_SUBMIT_DIR'),args['prefix'])
  else:
      # We are running in interactive mode
      args['workdir'] = os.path.join(os.getcwd(),args['prefix'])

  # Create and cd workdir
  if not os.path.exists(args['workdir']):
        os.makedirs(args['workdir'])
  os.chdir(args['workdir'])

  # Save all parameters in json file
  jsonfile = open(os.path.join(args['workdir'],'00_INFO.json'),'w')      
  json.dump(args,jsonfile, sort_keys = True, indent = 4)
  jsonfile.close()    
  
  # build job startup command
  test16bin = os.path.join(os.getenv('VFDIR'),'ValidationTests','test16')
  print test16bin
  cmd = '%s %s -p %s '%(args['mpiexec'],test16bin,args['prefix'])
  for k in sorted(args.keys()):
    if k not in ['gceff','mpiexec','workdir','prefix','gamg','extraopts']:
      val = ('%s'%args[k]).strip('[ ]').replace(' ','')
      cmd += '-%s %s '%(k,val)
  if args['gamg']:
    cmd += '-U_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 -U_mg_levels_ksp_type chebyshev -U_mg_levels_pc_type sor -U_pc_gamg_agg_nsmooths 1 -U_pc_gamg_threshold 0 -U_pc_gamg_type agg -U_pc_type gamg '
    #cmd += '-V_mg_levels_ksp_chebyshev_estimate_eigenvalues 0,0.1,0,1.1 -V_mg_levels_ksp_type chebyshev -V_mg_levels_pc_type sor -V_pc_gamg_agg_nsmooths 1 -V_pc_gamg_threshold 0 -V_pc_gamg_type agg -V_pc_type gamg '
  cmd += args['extraopts']  
      
  print "Now running %s"%cmd
  os.system(cmd)
  return 0
    
import sys  
if __name__ == "__main__":
	sys.exit(main())

