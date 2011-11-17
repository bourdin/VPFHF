/*
   VFFlow_DarcySteadyState.h
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
   (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
*/
/*
 VFFlow.h
 VF - Fluid Flow
 */

#ifndef VFFLOW_H
#define VFFLOW_H
#include "PetscFixes.h"

extern PetscErrorCode FlowSolverInitialize(VFCtx *ctx);
extern PetscErrorCode GetFlowProp(FlowProp *flowprop, 	FlowUnit flowunit, ResProp resprop);
extern PetscErrorCode SETFlowBC(FLOWBC *BC, FlowCases flowcase);
extern PetscErrorCode VecApplyFlowBC(Vec RHS, FLOWBC *BC, ResProp resprop);
extern PetscErrorCode MatApplyFlowBC(Mat K, DA da, FLOWBC *BC);
extern PetscErrorCode FlowMatnVecAssemble(Mat K, Vec RHS, VFFields *fields, VFCtx *ctx);
extern PetscErrorCode FLow_Vecg(PetscReal *Kg_local, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, FlowProp flowpropty);
extern PetscErrorCode FLow_Vecf(PetscReal *Kf_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c, FlowProp flowpropty);
extern PetscErrorCode FLow_MatD(PetscReal *Kd_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, FlowProp flowpropty);
extern PetscErrorCode FLow_MatB(PetscReal *KB_ele, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei, PetscInt c);
extern PetscErrorCode FLow_MatA(PetscReal *A_local, CartFE_Element3D *e,  PetscInt ek, PetscInt ej, PetscInt ei);
extern PetscErrorCode FlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VecApplyWellFlowBC(Vec RHS, VFCtx *ctx);
extern PetscErrorCode VFFlow_DarcySteadyState(VFCtx *ctx, VFFields *fields);
#endif /* VFFLOW_H */
