/*
 VFFlow_TSFEM.h
 Finite element implementation of  Darcy flow 
 (c) 2010-2013 K. Yoshioka
*/

#ifndef VFFLOW_TSFEM_H
#define VFFLOW_TSFEM_H

extern PetscErrorCode FEMTSFlowSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode FEMTSFlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FormIFunction_P(TS ts,PetscReal t,Vec VelnPress,Vec VelnPressdot,Vec Func,void *user);
extern PetscErrorCode FormIJacobian_P(TS ts,PetscReal t,Vec VelnPress,Vec VelnPressdot,PetscReal shift,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user);
extern PetscErrorCode FormFunction_P(TS ts,PetscReal t,Vec vec1,Vec Func,void *user);
extern PetscErrorCode FEMTSMonitor(TS ts,PetscInt timestep,PetscReal timevalue,Vec pressure,void*);
extern PetscErrorCode MatApplyTSPressureBC(Mat K,Mat Klhs,VFBC *bcP);
extern PetscErrorCode FlowFEMTSSolve(VFCtx *ctx,VFFields *fields);
#endif /* VFFLOW_TSFEM_H */
