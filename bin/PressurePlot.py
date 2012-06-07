#!/usr/bin/env python
import sys

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Plot energy evolution for test15.')
    parser.add_argument('inputfile',type=argparse.FileType('r'),nargs='*',help='Input file',default=sys.stdin)
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    return parser.parse_args()

def main():
    import pymef90
    import matplotlib
    import numpy as np
    import matplotlib.pyplot as plt
    options = parse()

    if options.outputfile != None:
      matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    for f in options.inputfile:
        print f.name
        p = np.loadtxt(f)    
        plt.plot(p[:,1],p[:,2],label=f.name)
        plt.grid()
        plt.legend(loc=0)
        plt.xlabel('V')
        plt.ylabel('p')
        plt.title('Pressure vs. injected volume')

    if options.outputfile != None:
      plt.savefig(options.outputfile)
    else:
      plt.show()


if __name__ == "__main__":
        sys.exit(main())