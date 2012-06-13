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
    #parser.add_argument('-o','--outputfile',help='output file',default=sys.stdout)
    parser.add_argument('-o','--outputfile', type=argparse.FileType('w', 0),default='-')
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
    
    #f = open(options.outputfile,'w')
    options.outputfile.write('PBS_JOBID')
    for k in sorted(allKeys):
        options.outputfile.write(',  ')
        options.outputfile.write(k)
    options.outputfile.write('\n')
    for D in allDict:
        options.outputfile.write(repr(D['PBS_JOBID']))
        for k in sorted(allKeys):
            options.outputfile.write(', ')
            try:
                options.outputfile.write(repr(D[k]))
            except KeyError:
                pass
        options.outputfile.write('\n')   

    #  plt.savefig(options.outputfile)
    #else:
    #  plt.show()


if __name__ == "__main__":
        sys.exit(main())
