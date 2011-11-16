/*
   VFFlow_DarcySteadyState.h
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
   (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
*/

#ifndef VFFLOW_MIXEDFEM_H
#define VFFLOW_MIXEDFEM_H
#include "PetscFixes.h"


extern PetscErrorCode VFFlow_MixedFEM(VFCtx *ctx, VFFields *fields);

#endif /* VFFLOW_MIXEDFEM_H */
