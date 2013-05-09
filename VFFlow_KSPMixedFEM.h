/*
   VFFlow_DarcySteadyState.h
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
 (c) 2011-2012 C. Chukwudozie, LSU
*/

#ifndef VFFLOW_KSPMIXEDFEM_H
#define VFFLOW_KSPMIXEDFEM_H

extern PetscErrorCode VFFlow_DarcyMixedFEMSteadyState(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode MixedFEMFlowSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode GetFlowProp(VFFlowProp *flowprop, VFUnit flowunit, VFResProp resprop);
extern PetscErrorCode ResetFlowBC(BC *bcP,BC *bcQ, VFFlowCases flowcase);
extern PetscErrorCode VecApplyVelocityBC(Vec RHS,BC *bcQ,VFCtx *ctx, PetscReal ****UnPre_array);
extern PetscErrorCode FlowMatnVecAssemble(Mat K, Mat Krhs, Vec RHS, VFFields *fields, VFCtx *ctx);
extern PetscErrorCode FLow_Vecg(PetscReal *Kg_local, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, FlowProp flowpropty, PetscReal ****perm_array, PetscReal ***v_array);
extern PetscErrorCode FLow_Vecf(PetscReal *Kf_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c, FlowProp flowpropty, PetscReal ****perm_array, PetscReal ***v_array);
extern PetscErrorCode FLow_MatD(PetscReal *Kd_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, FlowProp flowpropty, PetscReal ****perm_array, PetscReal ***v_array);
extern PetscErrorCode FLow_MatB(PetscReal *KB_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c, PetscReal ***v_array);
extern PetscErrorCode FLow_MatA(PetscReal *A_local, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscReal ***v_array);
extern PetscErrorCode MixedFEMFlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode BoundaryFlowRate(PetscReal *vel, PetscInt c, PetscInt i, PetscInt j, PetscInt k, PetscReal hi, PetscReal hj, PetscReal hk, FlowProp flowpropty);
extern PetscErrorCode BoundaryPressure(PetscReal *press, PetscInt i, PetscInt j, PetscInt k, PetscReal hi, PetscReal hj, PetscReal hk, ResProp resprop);
extern PetscErrorCode FLow_MatBTranspose(PetscReal *KB_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c, FlowProp flowpropty, PetscReal ****perm_array, PetscReal ***v_array);
extern PetscErrorCode VecApplySourceTerms(PetscReal *Ks_local, PetscReal ***source_array, CartFE_Element3D *e, PetscInt ek, PetscInt ej, PetscInt ei, VFCtx *ctx, PetscReal ***v_array);
extern PetscErrorCode ReSETSourceTerms(Vec Src, FlowProp flowpropty);
extern PetscErrorCode ReSETBoundaryTerms(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode MixedFEMKSPMonitor(KSP ksp,PetscInt its,PetscReal fnorm,void* ptr);
extern PetscErrorCode MatApplyKSPVelocityBC(Mat K,Mat Klhs,BC *bcQ);
extern PetscErrorCode VecApplyPressureBC(PetscReal *RHS_local,PetscReal ***pre_array,PetscInt ek,PetscInt ej,PetscInt ei,FACE face,CartFE_Element2D *e,FlowProp flowpropty,PetscReal ****perm_array, PetscReal ***v_array);
extern PetscErrorCode VecApplyWEllFlowRate(PetscReal *RHS_local,CartFE_Element3D *e,PetscReal Q,PetscReal hwx,PetscReal hwy,PetscReal hwz,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ***v_array);


#endif /* VFFLOW_MIXEDFEM_H */
