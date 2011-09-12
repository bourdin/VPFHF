/*
   VFFlow.h
   VF - Fluid Flow
*/

#ifndef VFU_H
#define VFU_H
#include "PetscFixes.h"

extern PetscErrorCode VFFlowTimeStep(VFCtx *ctx,VFFields *fields);

extern PetscErrorCode VF_FakeFlow(VFCtx *ctx, VFFields *fields,ResProp *resprop)


#endif /* VFFLOW_H */
