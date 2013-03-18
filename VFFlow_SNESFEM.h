/*
   VFFlow_SNESFEM.h
    
*/

#ifndef VFFLOW_SNESFEM_H
#define VFFLOW_SNESFEM_H

extern PetscErrorCode FormSNESMatricesnVector_P(Mat Kneu,Mat Kalt,Vec RHS,VFCtx *ctx);
extern PetscErrorCode FormSNESIFunction_P(SNES snes,Vec pressure,Vec Func,void *user);
extern PetscErrorCode FormSNESIJacobian_P(SNES snes,Vec pressure,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user);
extern PetscErrorCode FlowFEMSNESSolve(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FEMSNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void* ptr);
extern PetscErrorCode FEMSNESFlowSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode FEMSNESFlowSolverFinalize(VFCtx *ctx,VFFields *fields);

#endif /* VFFLOW_KSPFEM_H */
