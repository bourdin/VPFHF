/*
   VFFlow.h
   Generic interface to flow solvers

     (c) 2010-2012 Blaise Bourdin, LSU. bourdin@lsu.edu
                    Keita Yoshioka, Chevron ETC
*/

#if !defined(VFFLOW_H)
#define VFFLOW_H


typedef struct {
  PetscReal   mu;       /* Fluid viscosity          */
  PetscReal   rho;      /* Fluid density            */
  PetscReal   cf;       /* Fluid compressibility    */
  PetscReal   betac;    /* Conversion constant      */
  PetscReal   gammac;   /*Conversion parameter      */
  PetscReal   alphac;   /*Conversion parameter      */
  PetscReal   g;
} FluidProp;

extern PetscErrorCode FlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FlowSolverInitialize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode BCQInit(VFBC *BCQ,VFCtx *ctx);
extern PetscErrorCode BCFracQInit(VFBC *BCFracQ,VFCtx *ctx);
extern PetscErrorCode BCPInit(VFBC *BCP,VFCtx *ctx);
extern PetscErrorCode SETBoundaryTerms_P(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VF_StepP(VFFields *fields,VFCtx *ctx);
extern PetscErrorCode VecApplyPressureBC_FEM(Vec RHS,Vec BCF,VFBC *BC);
extern PetscErrorCode MatApplyPressureBC_FEM(Mat K,Mat M,VFBC *bcP);
extern PetscErrorCode VFFlow_FEM_MatKPAssembly3D_local(PetscReal *Mat_local,VFFlowProp *flowprop,PetscReal ****perm_array, PetscInt ek,PetscInt ej,PetscInt ei,VFCartFEElement3D *e);
extern PetscErrorCode VFFlow_FEM_MatMPAssembly3D_local(PetscReal *Mat_local,VFFlowProp *flowprop,PetscInt ek,PetscInt ej,PetscInt ei,PetscReal ACoef_P,VFCartFEElement3D *e);
extern PetscErrorCode Flow_Vecg_FEM(PetscReal *Kg_local,VFCartFEElement3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFFlowProp flowpropty,PetscReal ****perm_array);
extern PetscErrorCode VecApplySourceTerms_FEM(PetscReal *Ks_local,PetscReal ***source_array,VFCartFEElement3D *e,PetscInt ek,PetscInt ej,PetscInt ei,VFCtx *ctx);
extern PetscErrorCode VecApplyWellFlowRate_FEM(PetscReal *RHS_local,PetscReal Q,PetscReal hwx,PetscReal hwy,PetscReal hwz);
extern PetscErrorCode FEMSNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void* ptr);
extern PetscErrorCode GetFlowProp(VFFlowProp *flowprop,VFResProp *resprop,VFMatProp *matprop,VFCtx *ctx,VFFields *fields,PetscInt n);
extern PetscErrorCode VecApplyPressureBC_SNES(Vec Func,Vec pressure, Vec BCF,VFBC *BC);
#endif /* VFFLOW_H */

