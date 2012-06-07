#!/usr/bin/env python
import sys

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Plot energy evolution for VarFracQS.')
    parser.add_argument('inputfile',type=argparse.FileType('r'),nargs='?',help='Input file',default=sys.stdin)
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    parser.add_argument("-d","--debug",action="store_true",default=False,help="Display useless debugging information")
    parser.add_argument("-m","--stepmin",type=int,help="first time step")
    parser.add_argument("-M","--stepmax",type=int,help="last time step")
    return parser.parse_args()

def plot(energies):
  import matplotlib.pyplot as plt
  plt.plot(energies[:,0], energies[:,1], 'r-', label='Elastic energy')
  plt.plot(energies[:,0], energies[:,2], 'g-', label='Insitu work')
  plt.plot(energies[:,0], energies[:,3], 'b-', label='Surface energy')
  plt.plot(energies[:,0], energies[:,4], 'y-', label='Pressure work')
  plt.plot(energies[:,0], energies[:,1]+energies[:,2]+energies[:,3]+energies[:,4], 'k-', label='Total energy', lw=2)
  plt.grid()
  plt.legend(loc=0)
  plt.xlabel('Step')
  plt.ylabel('Energy')
  plt.title('Energies vs step')
  ###
  return 0


def main():
    import matplotlib
    import numpy as np
    options = parse()

    if options.debug:
      print("Option table: {0}".format(options))
    
    energies=np.loadtxt(options.inputfile)
      
    if options.stepmin == None:
      tmin = 0
    else:
      tmin = int(options.stepmin)
    if options.stepmax == None:
      tmax = energies.shape[0]
    else:
      tmax = int(options.stepmax)
    
    if options.debug:
      print('size of energies: {0}'.format(energies.shape))
      print('requested slice: {0}:{1}'.format(tmin,tmax))
      print('Energies: {0}'.format(energies[tmin:tmax,:]))
    
    if options.outputfile != None:
      matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    ### plot
    plot(energies[tmin:tmax,:])
      
    ### export plot if needed
    if options.outputfile != None:
      plt.savefig(options.outputfile)
    else:
      plt.show()

if __name__ == "__main__":
        sys.exit(main())
