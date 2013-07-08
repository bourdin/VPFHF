/*
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
 (c) 2011-2012 C. Chukwudozie, LSU
*/

#ifndef VFFLOW_FractureFlow_H
#define VFFLOW_FractureFlow_H

extern PetscErrorCode MixedFractureFlowSolverInitialize(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode MixedFractureFlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode MixedFracFlowSNESSolve(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FormFracSNESIJacobian(SNES snes,Vec VelnPress,Mat *Jac,Mat *Jacpre,MatStructure *str,void *user);
extern PetscErrorCode FormFracSNESIFunction(SNES snes,Vec VelnPress,Vec Func,void *user);
extern PetscErrorCode FormFracMatricesnVector(Mat K,Mat Klhs,Vec RHS,VFCtx *ctx, VFFields *fields);
extern PetscErrorCode FracFlow_MatB(PetscReal *KB_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt c,PetscReal ****u_array,PetscReal ***v_array);
extern PetscErrorCode FracFlow_MatBTranspose(PetscReal *KB_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscInt c,VFFlowProp flowpropty,PetscReal ****u_array,PetscReal ***v_array);
extern PetscErrorCode FracFlow_MatA(PetscReal *A_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ****u_array,PetscReal ***v_array);
extern PetscErrorCode FracFlow_Vecg(PetscReal *Kg_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ****u_array,PetscReal ***v_array);
extern PetscErrorCode FracFlow_MatS(PetscReal *S_local,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ****u_array,PetscReal ***v_array);
extern PetscErrorCode FracFlow_MatD(PetscReal *Kd_ele,CartFE_Element3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFFlowProp flowpropty,PetscReal ****u_array,PetscReal ***v_array);


#endif /* VFFLOW_MIXEDFEM_H */
