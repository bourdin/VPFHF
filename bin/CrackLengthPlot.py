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
    parser.add_argument('-d','--dim',help='Dimension',type=int,default=2)
    parser.add_argument('-t','--title',help='Title',default=None)
    return parser.parse_args()

def main():
    import matplotlib
    import numpy as np
    import os.path
    options = parse()

    if options.outputfile != None:
      matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    for f in options.inputfile:
        if (os.path.isfile(f)):
            p = np.loadtxt(f)   
            pp = p[:,3]
            if options.dim ==2:
                pp -= p[0,3]
            plt.plot(p[:,1],pp,'--',lw=2,label=f)
        if (os.path.isdir(f)):
            infotxt = os.path.join(f,'00_INFO.txt')
            if (os.path.exists(infotxt)):
                print 'Reading parameters for %s'%infotxt
                D = Dictreadtxt(infotxt)
                presfile = os.path.join(f,f+'.pres')
                if os.path.exists(presfile):
                    p = np.loadtxt(presfile)    
                    hx = D['LX']/(D['NX']+0.0)
                    l = '$h=%.2E$, $\epsilon/h=%.2f\ (%s)$'%(hx,D['EPSILON']/hx,f)
                    l = '$h=%.2E$, $\epsilon/h=%.2f$'%(hx,D['EPSILON']/hx)
                    #plt.plot(p[:,1]/D['LZ'],p[:,3]/D['LZ'],label=l,lw=2)
                    pp = p[:,3]
                    if options.dim ==2:
                        pp -= p[0,3]
                    plt.plot(p[:,1],pp,label=l,lw=2)
                
    plt.legend(loc=0)
    if options.dim ==2:
        plt.xlabel('Injected volume $(\mathrm{m}^2)$')
        plt.ylabel('Crack increment $(\mathrm{m})$')
    else:
        plt.xlabel('Injected volume $(\mathrm{m}^3)$')
        plt.ylabel('Crack radius $(\mathrm{m})$')

    ax = plt.gca()
    ax.grid()
    #ax.axis([0,.025,0,.015])
    plt.legend(loc=0,labelspacing=.1)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize='small')

    if not options.title:
        if options.dim ==2:
            plt.title('Crack increment vs. injected volume')
        else:
            plt.title('Crack radius vs. injected volume')
    else:
        plt.title(options.title)
        

    if options.outputfile != None:
      plt.savefig(options.outputfile)
    else:
      plt.show()


if __name__ == "__main__":
        sys.exit(main())
