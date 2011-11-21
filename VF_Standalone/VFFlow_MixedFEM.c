/*
   VFFlow_DarcySteadyState.c
   A mixed finite elements Darcy solver based on the method presented in
    [Chukwudi, please add reference here]
    
   (c) 2011 B. Bourdin.C. Chukwudozie, LSU, K. Yoshioka, CHEVRON ETC
*/

#include "petsc.h"
#include "CartFE.h"
#include "VFCommon.h"

#undef __FUNCT__
#define __FUNCT__ "VFFlow_MixFEM"
/* 
   Fake flow solver for VF_Chevron.c test
*/
extern PetscErrorCode VFFlow_MixedFEM(VFCtx *ctx, VFFields *fields)
{
  SETERRQ1(PETSC_ERR_SUP,"Flow solver %s not implemented yet",VFFlowSolverName[ctx->flowsolver]);
  PetscFunctionReturn(0);
}

