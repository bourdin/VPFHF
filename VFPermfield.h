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
extern PetscErrorCode CrackOpeningDisplacement(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode ComputeXYZOpening(CartFE_Element3D *e, PetscInt ei, PetscInt ej, PetscInt ek, PetscReal hx, PetscReal hy, PetscReal hz, PetscReal ****displ_array, PetscReal ***vfield_array, PetscReal ****perm_array);
extern PetscErrorCode NodeToCellInterpolation(DM dm, Vec node_vec, Vec cell_vec);

#endif 
