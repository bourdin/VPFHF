/*
 VFHeat.h
 Generic interface to heat solvers
 
 (c) 2011-2013 C. Chukwudozie, LSU
 */

#ifndef VFHEAT_H
#define VFHEAT_H

extern PetscErrorCode BCTInit(VFBC *BCT,VFCtx *ctx);
extern PetscErrorCode BCQTInit(VFBC *BCQT,VFCtx *ctx);
extern PetscErrorCode VF_HeatSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VF_HeatSolverInitialize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VF_HeatTimeStep(VFCtx *ctx,VFFields *fields);



#endif /* VFFLOW_H */

