/*
   VFFlow_DarcyPoisson.h
   A direct solver for the Darcy equation 
    
   (c) 2011 K. Yoshioka, CHEVRON ETC
*/

#ifndef VFFLOW_POISSON_H
#define VFFLOW_POISSON_H
#include "PetscFixes.h"


extern PetscErrorCode VFFlow_DarcyPoisson(VFCtx *ctx, VFFields *fields);

#endif /* VFFLOW_POISSON_H */
