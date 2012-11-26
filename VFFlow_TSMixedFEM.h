/*
   VFFlow_DarcySteadyState.h
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
 (c) 2011-2012 C. Chukwudozie, LSU
*/

#ifndef VFFLOW_TSMIXEDFEM_H
#define VFFLOW_TSMIXEDFEM_H

extern PetscErrorCode MixedFlowFEMTSSolve(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode MixedFEMTSFlowSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode MixedFEMTSFlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FormIFunction(TS ts,PetscReal t,Vec VelnPress,Vec VelnPressdot,Vec Func,void *user);
extern PetscErrorCode FormIJacobian(TS ts,PetscReal t,Vec VelnPress,Vec VelnPressdot,PetscReal shift,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user);
extern PetscErrorCode MatApplyTSVelocityBC(Mat K, Mat Klhs,FLOWBC *BC);
extern PetscErrorCode FormTSMatricesnVector(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx);
extern PetscErrorCode FormInitialSolution(Vec VelnPress,Vec VelnPressBV,FLOWBC *BC,VFCtx *ctx);
extern PetscErrorCode MixedFEMTSMonitor(TS ts,PetscInt timestep,PetscReal timevalue,Vec VelnPress,void*);
extern PetscErrorCode VecApplyTSVelocityBC(Vec RHS,Vec BCV, FLOWBC *BC,VFCtx *ctx);
extern PetscErrorCode FormFunction(TS ts,PetscReal t,Vec vec1,Vec Func,void *user);

#endif /* VFFLOW_MIXEDFEM_H */
