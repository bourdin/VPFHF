import fnmatch
import os

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
        plotpngscript='visit -cli -nowin -s ' + os.path.join(scriptpath,'PNGplot.py')

        jobs = data['jobs']
        print 'Project: {0}'.format(data['project'])
        #print args
        #print jobs
    
        for jobid in sorted(jobs.keys()):
            job = jobs[jobid]
            print '{0:_^40}'.format(jobid)
            prefix = os.path.join(job['path'],job['info']['JOBID'])
            try:
                prefix = os.path.join(job['path'],job['info']['JOBID'])
            except KeyError:
                print '00_INFO.txt file has no JOBID entry, skipping' 

            if args.plotener or args.all: 
                infile = prefix+'.pres'
                outfile = prefix+'_ener.pdf'
                if not check_mtime(infile, outfile):
                    print "Generating energy plot %s"%outfile
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
                    if (job['info']['NX'] == '2') or (job['info']['NY'] == '2') or (job['info']['NZ'] == '2'):
                        dim = 2
                    else:
                        dim = 3
                    cmd = 'cd %s; %s'%(job['path'],plotpngscript)
                    #print cmd
                    os.system(cmd)

