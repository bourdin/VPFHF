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
extern PetscErrorCode VF_FormFlowStandardFEMIJacobian(SNES snes,Vec Pressure,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user);

#endif 
