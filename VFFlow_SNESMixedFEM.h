/*
   VFFlow_DarcySteadyState.h
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
 (c) 2011-2012 C. Chukwudozie, LSU
*/

#ifndef VFFLOW_SNESMIXEDFEM_H
#define VFFLOW_SNESMIXEDFEM_H

extern PetscErrorCode MixedFlowFEMSNESSolve(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FormSNESMatricesnVector(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx);
extern PetscErrorCode FormSNESIJacobian(SNES snes,Vec VelnPress,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user);
extern PetscErrorCode FormSNESIFunction(SNES snes,Vec VelnPress,Vec Func,void *user);
extern PetscErrorCode MixedFEMSNESFlowSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode MixedFEMSNESFlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode MatApplySNESVelocityBC(Mat K,Mat Klhs,FLOWBC *BC);
extern PetscErrorCode MixedFEMSNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void* ptr);
extern PetscErrorCode VecApplySNESVelocityBC(Vec RHS,Vec BCV, FLOWBC *BC,VFCtx *ctx);
#endif /* VFFLOW_MIXEDFEM_H */
