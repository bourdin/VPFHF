/*
   VFFlow_DarcySteadyState.h
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
   (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
*/

#ifndef VFFLOW_MIXEDFEM_H
#define VFFLOW_MIXEDFEM_H

extern PetscErrorCode VFFlow_DarcyMixedFEMSteadyState(VFCtx *ctx, VFFields *fields);

/* 
  Rename and check if all these need to be public
*/
extern PetscErrorCode MixedFEMFlowSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode GetFlowProp(FlowProp *flowprop, 	FlowUnit flowunit, ResProp resprop);
extern PetscErrorCode SETFlowBC(FLOWBC *BC, FlowCases flowcase);
//extern PetscErrorCode VecApplyFlowBC(Vec RHS, FLOWBC *BC, VFCtx *ctx);
extern PetscErrorCode VecApplyFlowBC(Vec RHS,FLOWBC *BC,VFCtx *ctx, PetscReal ****UnPre_array);
extern PetscErrorCode MatApplyFlowBC(Mat K, DM da, FLOWBC *BC);
extern PetscErrorCode FlowMatnVecAssemble(Mat K, Mat Krhs, Vec RHS, VFFields *fields, VFCtx *ctx);
extern PetscErrorCode FLow_Vecg(PetscReal *Kg_local, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, FlowProp flowpropty, PetscReal ****perm_array);
extern PetscErrorCode FLow_Vecf(PetscReal *Kf_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c, FlowProp flowpropty, PetscReal ****perm_array);
extern PetscErrorCode FLow_MatD(PetscReal *Kd_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, FlowProp flowpropty, PetscReal ****perm_array);
extern PetscErrorCode FLow_MatB(PetscReal *KB_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c);
extern PetscErrorCode FLow_MatA(PetscReal *A_local, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei);
extern PetscErrorCode MixedFEMFlowSolverFinalize(VFCtx *ctx,VFFields *fields);
//extern PetscErrorCode VecApplyWellFlowBC(Vec RHS, VFCtx *ctx);
extern PetscErrorCode BoundaryFlowRate(PetscReal *vel, PetscInt c, PetscInt i, PetscInt j, PetscInt k, PetscReal hi, PetscReal hj, PetscReal hk, FlowProp flowpropty);
extern PetscErrorCode BoundaryPressure(PetscReal *press, PetscInt i, PetscInt j, PetscInt k, PetscReal hi, PetscReal hj, PetscReal hk, ResProp resprop);
extern PetscErrorCode FLow_MatBTranspose(PetscReal *KB_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c, FlowProp flowpropty, PetscReal ****perm_array);
extern PetscErrorCode VecApplyWellFlowBC(PetscReal *Ks_local, PetscReal ***source_array, CartFE_Element3D *e, PetscInt ek, PetscInt ej, PetscInt ei, VFCtx *ctx);
extern PetscErrorCode SETSourceTerms(Vec Src, FlowProp flowpropty);
extern PetscErrorCode SETBoundaryTerms(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode MixedFEMKSPMonitor(KSP ksp,PetscInt its,PetscReal fnorm,void* ptr);
extern PetscErrorCode ApplyKSPJacobianBC(Mat K,Mat Klhs,FLOWBC *BC);

/*	TS Functions*/
extern PetscErrorCode MixedFEMTSFlowSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode MixedFEMTSFlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FormIFunction(TS ts,PetscReal t,Vec VelnPress,Vec VelnPressdot,Vec Func,void *user);
extern PetscErrorCode FormIJacobian(TS ts,PetscReal t,Vec VelnPress,Vec VelnPressdot,PetscReal shift,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user);
extern PetscErrorCode ApplyTSJacobianBC(Mat K, Mat Klhs,FLOWBC *BC);
extern PetscErrorCode FormTSMatricesnVector(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx);
extern PetscErrorCode MixedFlowFEMTSSolve(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FormInitialSolution(Vec VelnPress,Vec VelnPressBV,FLOWBC *BC,VFCtx *ctx);
extern PetscErrorCode MixedFEMTSMonitor(TS ts,PetscInt timestep,PetscReal timevalue,Vec VelnPress,void*);
extern PetscErrorCode VecApplyTSFlowBC(Vec RHS,Vec bcv_local, FLOWBC *BC,VFCtx *ctx);
/*	SNES Functions*/
extern PetscErrorCode FormSNESMatricesnVector(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx);
extern PetscErrorCode FormSNESIJacobian(SNES snes,Vec VelnPress,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user);
extern PetscErrorCode FormSNESIFunction(SNES snes,Vec VelnPress,Vec Func,void *user);
extern PetscErrorCode MixedFlowFEMSNESSolve(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode MixedFEMSNESFlowSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode MixedFEMSNESFlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode ApplySNESJacobianBC(Mat K, Mat Klhs,FLOWBC *BC);
extern PetscErrorCode MixedFEMSNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void* ptr);









#endif /* VFFLOW_MIXEDFEM_H */
