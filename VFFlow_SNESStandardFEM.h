/*
 VFFlow_SNESStandardFEM.h
 A standard finite elements Darcy solver 
 [Chukwudi, please add reference here]
 
 (c) 2011-2013 C. Chukwudozie, LSU
 */

#ifndef VFFLOW_SNESSTANDARDDFEM_H
#define VFFLOW_SNESSTANDARDDFEM_H

extern PetscErrorCode VFFlow_SNESStandardFEMInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VFFlow_SNESStandardFEMFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VF_FormFlowStandardFEMIFunction(SNES snes,Vec Pressure,Vec Func,void *user);
extern PetscErrorCode VF_FlowStandardFEMSNESSolve(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VF_FormFlowStandardFEMMatricesnVectors(Mat K,Mat Krhs,Vec RHS,VFFields * fields,VFCtx *ctx);
extern PetscErrorCode VF_FormFlowStandardFEMIJacobian(SNES snes,Vec Pressure,Mat Jac,Mat Jacpre,void *user);
extern PetscErrorCode VF_FlowRateCompute(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VF_FlowRateCompute_local(PetscReal ****cellflowrate_array, PetscReal ***press_array ,PetscReal ****perm_array, PetscReal ***v_array, VFFlowProp *flowpropty, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
extern PetscErrorCode FlowVelocityCompute(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode FlowVelocityCompute_local(PetscReal ****cellvelocityrate_array, PetscReal ***press_array, PetscReal ****perm_array, PetscReal ***v_array, VFFlowProp *flowpropty, PetscInt ek, PetscInt ej, PetscInt ei, VFCartFEElement3D *e);
#endif 
