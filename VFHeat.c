/*
 VFHeat.c
 Generic interface to heat solvers
 
 (c) 2011-2013 C. Chukwudozie, LSU
 
 */
#include "petsc.h"
#include "VFCartFE.h"
#include "VFCommon.h"
#include "VFHeat.h"
#include "VFHeat_SNESFEM.h"



#undef __FUNCT__
#define __FUNCT__ "BCTInit"
extern PetscErrorCode BCTInit(BC *BCT,VFCtx *ctx)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = BCInit(BCT,1);CHKERRQ(ierr);	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BCQTInit"
extern PetscErrorCode BCQTInit(BC *BCQT,VFCtx *ctx)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = BCInit(BCQT,3);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}


/*
 VFHeatTimeStep: Does one time step of the flow solver selected in ctx.heatsolver
 */
#undef __FUNCT__
#define __FUNCT__ "VF_HeatTimeStep"
extern PetscErrorCode VF_HeatTimeStep(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode     ierr;
	
	PetscFunctionBegin;
	switch (ctx->heatsolver) {
		case HEATSOLVER_SNESFEM:
			ierr = VF_HeatFEMSNESSolve(ctx,fields);CHKERRQ(ierr);
			break;
		default:
		  break;
	}
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VF_HeatSolverFinalize"
extern PetscErrorCode VF_HeatSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	switch (ctx->heatsolver) {
		case HEATSOLVER_SNESFEM:       
			ierr = VF_FEMSNESHeatSolverFinalize(ctx,fields);CHKERRQ(ierr);
			break;
		default:
		  break;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VF_HeatSolverInitialize"
extern PetscErrorCode VF_HeatSolverInitialize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	switch (ctx->heatsolver) {
		case HEATSOLVER_SNESFEM:       
			ierr = VF_FEMSNESHeatSolverInitialize(ctx,fields);CHKERRQ(ierr);
			break;
		default:
		  break;
	}
	PetscFunctionReturn(0);
}


















