# VPFHF: Variational Phase-Field approach to Hydraulic Fracturing

An implementation of the algorithms described in 
   * Chukwudozie, C. (2016). Application of the Variational Fracture Model to Hydraulic Fracturing in Poroelastic Media. PhD thesis, Louisiana State University, Craft & Hawkins Department of Petroleum Engineering.
   * Chukwudozie, C., Bourdin, B., and Yoshioka, K. (2019). A variational phase-field model for hydraulic fracturing in porous media. Comp. Meth. Appl. Mech. Engng., 347:957–982.

## Installation:
This software requires PETSc 3.10 or later. The environment variable ```VFDIR``` must point to the root directory of the folder. Upon executing ```make```, a binary ```VPFHF``` should be buiult in ```$VFDIR/bin/$PETSC_ARCH```

## Testing
Some examples are included in ```$VFDIR/ValidationTests```