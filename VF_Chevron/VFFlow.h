/*
   VFFlow.h
   VF - Fluid Flow
*/

#ifndef VFFLOW_H
#define VFFLOW_H
#include "PetscFixes.h"

extern PetscErrorCode VFFlowTimeStep(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VF_FakeFlow(VFCtx *ctx, VFFields *fields);


#endif /* VFFLOW_H */
