/*
   VFFlow.h
   Generic interface to flow solvers

     (c) 2010-2012 Blaise Bourdin, LSU. bourdin@lsu.edu
                    Keita Yoshioka, Chevron ETC
*/

#if !defined(VFFLOW_H)
#define VFFLOW_H


typedef struct {
  PetscReal       mu;   /* Fluid viscosity          */
  PetscReal       rho;  /* Fluid density            */
  PetscReal       cf;   /* Fluid compressibility    */
  PetscReal   betac;    /* Conversion constant      */
  PetscReal   gammac;   /*Conversion parameter      */
  PetscReal   alphac;   /*Conversion parameter      */
  PetscReal   g;
} FluidProp;

extern PetscErrorCode FlowSolverFinalize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode FlowSolverInitialize(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode BCQInit(BC *BCQ,VFCtx *ctx);
extern PetscErrorCode BCFracQInit(BC *BCFracQ,VFCtx *ctx);
extern PetscErrorCode BCPInit(BC *BCP,VFCtx *ctx);
extern PetscErrorCode SETBoundaryTerms_P(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode VFFlowTimeStep(VFCtx *ctx,VFFields *fields);
extern PetscErrorCode VecApplyPressureBC_FEM(Vec RHS,Vec BCF,BC *BC);
extern PetscErrorCode MatApplyPressureBC_FEM(Mat K,Mat M,BC *bcP);
extern PetscErrorCode VFFlow_FEM_MatKPAssembly3D_local(PetscReal *Mat_local,VFFlowProp *flowprop,PetscReal ****perm_array, PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e);
extern PetscErrorCode VFFlow_FEM_MatMPAssembly3D_local(PetscReal *Mat_local,VFFlowProp *flowprop,PetscInt ek,PetscInt ej,PetscInt ei,CartFE_Element3D *e);
extern PetscErrorCode FEMSNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void* ptr);
extern PetscErrorCode GetFlowProp(VFFlowProp *flowprop,VFUnit flowunit,VFResProp resprop);
extern PetscErrorCode ResetFlowBC(BC *bcP,BC *bcQ, VFFlowCases flowcase);
extern PetscErrorCode ResetBoundaryTerms(VFCtx *ctx, VFFields *fields);
extern PetscErrorCode ResetSourceTerms(Vec Src,VFFlowProp flowpropty);
#endif /* VFFLOW_H */

