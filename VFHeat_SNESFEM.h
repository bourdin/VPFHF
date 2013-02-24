/*
 VFHeat_SNESFEM.h
 
 (c) 2011-2013 C. Chukwudozie, LSU
 */

#ifndef VFHEAT_SNESFEM_H
#define VFHEAT_SNESFEM_H

extern PetscErrorCode FEMSNESHeatSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode FEMSNESHeatSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode HeatFEMSNESSolve(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode HeatFEMSNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void* ptr);
extern PetscErrorCode FormHeatSNESIJacobian(SNES snes,Vec VelnPress,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user);


#endif /* VFHEAT_SNESFEM_H */
