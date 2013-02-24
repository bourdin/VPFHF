/*
   VFFlow.c
   Generic interface to heat solvers

 (c) 2011-2013 C. Chukwudozie, LSU

*/
#include "petsc.h"
#include "CartFE.h"
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

/*
 VFHeatTimeStep: Does one time step of the flow solver selected in ctx.heatsolver
 */
#undef __FUNCT__
#define __FUNCT__ "VFHeatTimeStep"
extern PetscErrorCode VFHeatTimeStep(VFCtx *ctx,VFFields *fields)
{
	char           filename[FILENAME_MAX];
	PetscViewer    viewer;
	PetscErrorCode     ierr;
	
	PetscFunctionBegin;
	switch (ctx->heatsolver) {
		case HEATSOLVER_SNESFEM:
			ierr = HeatFEMSNESSolve(ctx,fields);CHKERRQ(ierr);
			break;
	}
	PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "HeatSolverFinalize"
extern PetscErrorCode HeatSolverFinalize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	switch (ctx->heatsolver) {
		case HEATSOLVER_SNESFEM:       
			ierr = FEMSNESHeatSolverFinalize(ctx,fields);CHKERRQ(ierr);
			break;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "HeatSolverInitialize"
extern PetscErrorCode HeatSolverInitialize(VFCtx *ctx,VFFields *fields)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	switch (ctx->heatsolver) {
		case HEATSOLVER_SNESFEM:       
			ierr = FEMSNESHeatSolverInitialize(ctx,fields);CHKERRQ(ierr);
			break;
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BCqInit"
extern PetscErrorCode BCqInit(BC *BCq,VFCtx *ctx)
{
	PetscErrorCode ierr;
	
	PetscFunctionBegin;
	ierr = BCInit(BCq,3);CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

















