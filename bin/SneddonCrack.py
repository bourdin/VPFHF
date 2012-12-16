#!/usr/bin/env python

def SneddonPressure2D(V,l0,lz,Gc,E,nu):
    import math
    Ep = E/(1.-nu*nu)
    V0 = math.sqrt(4.*math.pi*Gc*l0*l0*l0/Ep)*lz
    if V <= V0:
        p = V*Ep/2./math.pi/l0/l0/lz
    else:
        p = math.pow(2.*Ep*Gc*Gc*lz/math.pi/V, 1./3.)
    return p

def SneddonLength2D(V,l0,lz,Gc,E,nu):
    import math
    Ep = E/(1.-nu*nu)
    V0 = math.sqrt(4.*math.pi*Gc*l0*l0*l0/Ep)*lz
    if V <= V0:
        l = l0
    else:
        l = pow(V*V*Ep/4./math.pi/Gc/lz/lz,1./3.)
    return l

def SneddonPressure3D(V,l0,Gc,E,nu):
    import math
    Ep = E/(1.-nu*nu)
    V0 = math.sqrt(64.*math.pi*Gc*math.pow(l0,5)/Ep/9)
    if V <= V0:
        p = 2
    else:
        p = math.pow(Ep*Ep*Gc*Gc*Gc*math.pi*math.pi*math.pi/12./V, 1./5.)
    return p

def SneddonLength3D(V,l0,Gc,E,nu):
    import math
    Ep = E/(1.-nu*nu)
    V0 = math.sqrt(64.*math.pi*Gc*math.pow(l0,5)/Ep/9)
    if V <= V0:
        l = l0
    else:
        l = pow(9.*V*V*Ep/64./math.pi/Gc,1./3.)
    return l

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Plot energy evolution for test15.')
    parser.add_argument('-o','--outputfile', type=argparse.FileType('w', 0),default='-')
    parser.add_argument('-e',type=float,help="Young's modulus",default=1.)
    parser.add_argument('-n','--nu',type=float,help="Poisson ratio",default=0.)
    parser.add_argument('--gc',type=float,help="fracture toughness",default=1.)
    parser.add_argument('--vmin',type=float,help="Min injected volume",default=0.)
    parser.add_argument('--vmax',type=float,help="Max injected volume",default=.03)
    parser.add_argument('--l0',type=float,help="initial crack length",default=.2)
    parser.add_argument('--lz',type=float,help="thickness",default=1)
    parser.add_argument('--dim',type=int,help="dimension",choices=[2,3],default=2)
    return parser.parse_args()

def main():
    import matplotlib
    import numpy as np
    import matplotlib.pyplot as plt
    import os.path
    options = parse()

    n = 151
    V = np.linspace(options.vmin,options.vmax,n)

    options.outputfile.write("# i   \tp \t l\n")
    for i in range(n):
        if options.dim == 2:
            p = SneddonPressure2D(V[i],options.l0,options.lz,options.gc,options.e,options.nu)
            l = options.gc*SneddonLength2D(V[i],options.l0,options.lz,options.gc,options.e,options.nu)
            options.outputfile.write('%i %e %e %e\n'%(i,V[i],p,2.*l*options.gc*options.lz))
        else:
            p = SneddonPressure3D(V[i],options.l0,options.gc,options.e,options.nu)
            l = options.gc*SneddonLength3D(V[i],options.l0,options.gc,options.e,options.nu)
            #options.outputfile.write('%i %e %e %e\n'%(i,V[i],p,l*options.gc))
            options.outputfile.write('%i %e %e %e\n'%(i,V[i],p,l))

    plt.show()

if __name__ == "__main__":
    import sys
    sys.exit(main())
