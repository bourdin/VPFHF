import fnmatch
import os

npmax = 5
from . import SNEDDONAction, scan, check_mtime
class all(SNEDDONAction):
    def action(self, args):
        import os.path
        data = scan(args.outdir,args.paths,args.excludedirs,
        args.excludefiles,args.project)
        
        scriptpath = os.path.dirname(__file__)
        energyplotscript = os.path.join(os.getenv('VFDIR'),'bin','plotener.py')
        lengthplotscript = os.path.join(os.getenv('VFDIR'),'bin','CrackLengthPlot.py')
        pressureplotscript = os.path.join(os.getenv('VFDIR'),'bin','PressurePlot2.py')
        pngscript = 'visit -cli -nowin -s ' + os.path.join(scriptpath,'PNGplot.py')
        pngscriptdist = os.path.join(os.path.dirname(__file__),'PNGplot.sh')
        movtransientscript = 'visit -cli -nowin -s ' + os.path.join(scriptpath,'MOVTransient.py')
        movtransientscriptdist = os.path.join(os.path.dirname(__file__),'MOVTransient.sh')
        mov3Dscript = 'visit -cli -nowin -s ' + os.path.join(scriptpath,'MOV3D.py')
        mov3Dscriptdist = os.path.join(os.path.dirname(__file__),'MOV3D.sh')
        jobs = data['jobs']
        print 'Project: {0}'.format(data['project'])


        for jobid in sorted(jobs.keys()):
            job = jobs[jobid]
            print '{0:_^40}'.format(jobid)
            prefix = os.path.join(job['path'],job['info']['JOBID'])
            try:
                prefix = os.path.join(job['path'],job['info']['JOBID'])
            except KeyError:
                print '00_INFO.txt file has no JOBID entry, skipping' 

            nx = float(job['info']['NX'])
            ny = float(job['info']['NY'])
            nz = float(job['info']['NZ'])
            ne = int(nx * ny * nz)
            N = min(max(1,ne / 300000 / 12),npmax) # limit to 400,000 elements / cpu
            n = N * 12

            if args.plotener: 
                infile = prefix+'.pres'
                outfile = prefix+'_ener.pdf'
                if not check_mtime(infile, outfile):
                    print "Generating energy plot   %s"%outfile
                    cmd = 'python %s %s -o %s'%(energyplotscript,infile,outfile)
                    os.system(cmd)

            if args.plotlength or args.all: 
                infile = prefix+'.pres'
                outfile = prefix+'_length.pdf'
                if not check_mtime(infile, outfile):
                    print "Generating length plot   %s"%outfile
                    if (job['info']['NX'] == '2') or (job['info']['NY'] == '2') or (job['info']['NZ'] == '2'):
                        dim = 2
                    else:
                        dim = 3
                    cmd = 'python %s %s -o %s --dim %i'%(lengthplotscript,infile,outfile,dim)
                    os.system(cmd)

            if args.plotpressure or args.all: 
                infile = prefix+'.pres'
                outfile = prefix+'_pres.pdf'
                if not check_mtime(infile, outfile):
                    print "Generating pressure plot %s"%outfile
                    if (job['info']['NX'] == '2') or (job['info']['NY'] == '2') or (job['info']['NZ'] == '2'):
                        dim = 2
                    else:
                        dim = 3
                    cmd = 'python %s %s -o %s --dim %i'%(pressureplotscript,infile,outfile,dim)
                    os.system(cmd)

            if args.plotpng or args.all: 
                infile = prefix+'.xmf'
                outfile = prefix+'.png'
                if not check_mtime(infile, outfile):
                    print "Generating png plot %s"%outfile
                    cmd = 'cd %s; %s'%(job['path'],pngscript)
                    #print cmd
                    os.system(cmd)

            if args.movtransient: 
                infile = prefix+'.xmf'
                outfile = prefix+'-Transient.avi'
                if not check_mtime(infile, outfile):
                    print "Generating transient movie frames"
                    if args.dist:
                        cmd = 'sbatch -N%i -n%i -J %s %s %s'%(N,n,job['info']['JOBID'],movtransientscriptdist,job['path'])
                    else:
                        cmd = 'cd %s; %s'%(job['path'],movtransientscript)
                    #print cmd
                    os.system(cmd)

            if args.mov3d: 
                infile = prefix+'.xmf'
                outfile = prefix+'-3D.avi'
                if not check_mtime(infile, outfile):
                    print "Generating 3D movie frames"
                    if args.dist:
                        cmd = 'sbatch -N%i -n%i -J %s %s %s'%(N,n,job['info']['JOBID'],mov3Dscriptdist,job['path'])
                    else:
                        cmd = 'cd %s; %s'%(job['path'],mov3Dscript)
                    os.system(cmd)

