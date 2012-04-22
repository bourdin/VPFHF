/*
   VFFlow_DarcySteadyState.h
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
   (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
*/

#ifndef VFPermfield_H
#define VFPermfield_H
#include "PetscFixes.h"


/* 
  Rename and check if all these need to be public
*/
extern PetscErrorCode VolumetricCrackOpening(PetscReal *CrackVolume, VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VolumetricCrackOpening3D_local(PetscReal ***volcrackopening_array, PetscReal ****displ_array, PetscReal ***vfield_array, PetscInt ek, PetscInt ej, PetscInt ei, CartFE_Element3D *e);
extern PetscErrorCode CellToNodeInterpolation(DM dm, Vec node_vec, Vec cell_vec, VFCtx *ctx);
#endif 
