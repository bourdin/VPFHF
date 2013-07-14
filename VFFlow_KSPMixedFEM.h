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
extern PetscErrorCode Flow_Vecg(PetscReal *Kg_local, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, VFFlowProp flowpropty, PetscReal ****perm_array, PetscReal ***v_array);
extern PetscErrorCode Flow_Vecf(PetscReal *Kf_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c, VFFlowProp flowpropty, PetscReal ***v_array);
extern PetscErrorCode Flow_MatD(PetscReal *Kd_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, VFFlowProp flowpropty, PetscReal ****perm_array, PetscReal ***v_array);
extern PetscErrorCode Flow_MatB(PetscReal *KB_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c, PetscReal ***v_array);
extern PetscErrorCode Flow_MatA(PetscReal *A_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei, PetscInt c,VFFlowProp flowpropty,PetscReal ****perm_array,PetscReal ***v_array);
extern PetscErrorCode MixedFEMFlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode BoundaryFlowRate(PetscReal *vel, PetscInt c, PetscInt i, PetscInt j, PetscInt k, PetscReal hi, PetscReal hj, PetscReal hk, VFFlowProp flowpropty);
extern PetscErrorCode BoundaryPressure(PetscReal *press, PetscInt i, PetscInt j, PetscInt k, PetscReal hi, PetscReal hj, PetscReal hk, VFResProp resprop);
extern PetscErrorCode Flow_MatBTranspose(PetscReal *KB_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c, PetscReal ***v_array);
extern PetscErrorCode VecApplySourceTerms(PetscReal *Ks_local, PetscReal ***source_array, CartFE_Element3D *e, PetscInt ek, PetscInt ej, PetscInt ei, VFCtx *ctx, PetscReal ***v_array);
extern PetscErrorCode ResetSourceTerms(Vec Src, VFFlowProp flowpropty);
extern PetscErrorCode ResetBoundaryTerms(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode MixedFEMKSPMonitor(KSP ksp,PetscInt its,PetscReal fnorm,void* ptr);
extern PetscErrorCode MatApplyKSPVelocityBC(Mat K,Mat Klhs,BC *bcQ);
extern PetscErrorCode VecApplyPressureBC(PetscReal *RHS_local,PetscReal ***pre_array,PetscInt ek,PetscInt ej,PetscInt ei,FACE face,CartFE_Element2D *e,VFFlowProp flowpropty,PetscReal ****perm_array, PetscReal ***v_array);
extern PetscErrorCode VecApplyWellFlowRate(PetscReal *RHS_local,CartFE_Element3D *e,PetscReal Q,PetscReal hwx,PetscReal hwy,PetscReal hwz,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ***v_array);
extern PetscErrorCode MixedFlowFEMKSPSolve(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VF_MatA_local(PetscReal *A_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ***v_array);
extern PetscErrorCode VF_RHSFlowMechUCoupling_local(PetscReal *K_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFFlowProp flowpropty,PetscReal ****u_diff_array,PetscReal ***v_array);
extern PetscErrorCode VF_RHSAddMechUCoupling(Vec RHS,VFCtx *ctx);

#endif /* VFFLOW_MIXEDFEM_H */
