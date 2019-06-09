# VPFHF: Variational Phase-Field approach to Hydraulic Fracturing

An implementation of the algorithms described in 
   * Chukwudozie, C. (2016). Application of the Variational Fracture Model to Hydraulic Fracturing in Poroelastic Media. PhD thesis, Louisiana State University, Craft & Hawkins Department of Petroleum Engineering.
   * Chukwudozie, C., Bourdin, B., and Yoshioka, K. (2019). A variational phase-field model for hydraulic fracturing in porous media. Comp. Meth. Appl. Mech. Engng., 347:957–982.

## License:
This software is free software, distributed under the 2-clause BSD license. A copy of the license is included in the LICENSE file.
We cordially ask that any published work derived from this code, or utilizing it references the software as [![DOI](https://zenodo.org/badge/191032060.svg)](https://zenodo.org/badge/latestdoi/191032060), and the above-mentioned published works:
```

@article{Chukwudozie-Bourdin-EtAl-2019a,
	Author = {Chukwudozie, C. and Bourdin, B. and Yoshioka, K.},
	Journal = {Comp. Meth. Appl. Mech. Engng.},
	Pages = {957--982},
	Title = {A variational phase-field model for hydraulic fracturing in porous media},
	Volume = {347},
	Year = {2019}}

@phdthesis{Chukwudozie-2016a,
	Address = {Craft \& Hawkins Department of Petroleum Engineering},
	Author = {Chukwudozie, C.},
	School = {Louisiana State University},
	Title = {Application of the Variational Fracture Model to Hydraulic Fracturing in Poroelastic Media},
	Year = {2016}}
```

## Installation:
This software requires PETSc 3.10 or later. The environment variable ```VFDIR``` must point to the root directory of the folder. Upon executing ```make```, a binary ```VPFHF``` should be buiult in ```$VFDIR/bin/$PETSC_ARCH```

## Testing:
Some examples are included in ```$VFDIR/ValidationTests```
