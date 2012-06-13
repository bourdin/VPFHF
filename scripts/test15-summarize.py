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
    parser = argparse.ArgumentParser(description='Creates a summary of a series of computations in csv format.')
    parser.add_argument('inputfile',nargs='*',help='Input file',default=sys.stdin)
    parser.add_argument('-o','--outputfile',help='output file',default=sys.stdout)
    return parser.parse_args()

def main():
    import os.path
    options = parse()

    allDict = []
    allKeys = set()
    for f in options.inputfile:
        if (os.path.isdir(f)):
            infotxt = os.path.join(f,'00_INFO.txt')
            if (os.path.exists(infotxt)):
                D = Dictreadtxt(infotxt)
                for k in D.keys():
                    allKeys.add(k)
                allDict.append(D)        
    allKeys.remove('PBS_JOBID')
    
    f = open(options.outputfile,'w')
    f.write('PBS_JOBID')
    for k in sorted(allKeys):
        f.write(',  ')
        f.write(k)
    f.write('\n')
    for D in allDict:
        f.write(repr(D['PBS_JOBID']))
        for k in sorted(allKeys):
            f.write(', ')
            try:
                f.write(repr(D[k]))
            except KeyError:
                pass
        f.write('\n')   

    #  plt.savefig(options.outputfile)
    #else:
    #  plt.show()


if __name__ == "__main__":
        sys.exit(main())
