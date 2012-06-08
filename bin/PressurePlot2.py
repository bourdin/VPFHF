#!/usr/bin/env python
import sys

def Dictreadtxt(filename):
    D = {}
    for l in open(filename,'r').readlines():
        l = l.strip()
        k, v = l.split(' ', 1) if l.count(' ') > 0 else (l, '')
        v = v.strip()
        try:
            v = float(v)            
            if int(v) == v:
                v = int(v)
        except ValueError:
            pass
        D[k] = v
    return D

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Plot energy evolution for test15.')
    #parser.add_argument('inputfile',type=argparse.FileType('r'),nargs='*',help='Input file',default=sys.stdin)
    parser.add_argument('inputfile',nargs='*',help='Input file',default=sys.stdin)
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    return parser.parse_args()

def main():
    import pymef90
    import matplotlib
    import numpy as np
    import matplotlib.pyplot as plt
    import os.path
    options = parse()
    print options

    if options.outputfile != None:
      matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    for f in options.inputfile:
        if (os.path.isfile(f)):
            print '%s is a file'%f
            p = np.loadtxt(f)    
            plt.plot(p[:,1],p[:,2],label=f)
        if (os.path.isdir(f)):
            print '%s is a folder'%f
            infotxt = os.path.join(f,'00_INFO.txt')
            if (os.path.exists(infotxt)):
                print 'Reading parameters for %s'%infotxt
                D = Dictreadtxt(infotxt)
                presfile = os.path.join(f,f+'.pres')
                if os.path.exists(presfile):
                    p = np.loadtxt(presfile)    
                    hx = D['LX']/(D['NX']+0.0)
                    l = '%s: $h=%.2E$, $\epsilon/h=%.1f$'%(f,hx,D['EPSILON']/hx)
                    plt.plot(p[:,1],p[:,2],lw=2,label=l)
                
    ax = plt.gca()
    ax.grid()
    ax.axis([0,.03,0,2.5])
    plt.legend(loc=0,labelspacing=.1)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize='small')

    plt.xlabel('V')
    plt.ylabel('p')
    plt.title('Pressure vs. injected volume')

    if options.outputfile != None:
      plt.savefig(options.outputfile)
    else:
      plt.show()


if __name__ == "__main__":
        sys.exit(main())
