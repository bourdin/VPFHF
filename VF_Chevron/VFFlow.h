/*
   VFFlow.h
   VF - Fluid Flow
*/

#ifndef VFFLOW_H
#define VFFLOW_H
#include "PetscFixes.h"

extern PetscErrorCode VFFlowTimeStep(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VF_FakeFlow(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VF_LinFlow(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VF_PAssembly3D(Mat K,Vec RHS,VFFields *fields,VFCtx *ctx);
extern PetscErrorCode VF_MatP3D_local(PetscReal *Mat_local,ResProp *resprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e);

#endif /* VFFLOW_H */
