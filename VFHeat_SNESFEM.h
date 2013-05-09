/*
 VFHeat_SNESFEM.h
 
 (c) 2011-2013 C. Chukwudozie, LSU
 */

#ifndef VFHEAT_SNESFEM_H
#define VFHEAT_SNESFEM_H

extern PetscErrorCode VF_FEMSNESHeatSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VF_FEMSNESHeatSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VF_HeatFEMSNESSolve(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FormSNESHeatIJacobian(SNES snes,Vec T,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user);
extern PetscErrorCode FormSNESHeatIFunction(SNES snes,Vec T,Vec Func,void *user);
extern PetscErrorCode FormHeatMatricesnVector(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx, VFFields *field);
extern PetscErrorCode VecApplyHeatFluxBC(PetscReal *RHS_local,PetscReal ***flux_array,PetscInt ek,PetscInt ej,PetscInt ei,FACE face,CartFE_Element2D *e,PetscReal ***v_array);
extern PetscErrorCode HeatMatN(PetscReal *KN_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei, PetscReal ****vel_array,PetscReal ***v_array);
extern PetscErrorCode VecApplyHeatSourceTerms(PetscReal *K1_local,PetscReal *K2_local,PetscReal ***source_array,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFCtx *ctx, PetscReal ****vel_array,PetscReal ***v_array);
extern PetscErrorCode HeatMatC(PetscReal *KC_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei, PetscReal ****vel_array,PetscReal ***v_array);



#endif /* VFHEAT_SNESFEM_H */
